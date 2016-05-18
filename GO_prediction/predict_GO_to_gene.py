"""
    File for predicting whether a given GO term is associated with
    various genes.

    Consider biological process X
    Training example: gene
    Label (binary): 1 if gene is associated with X in GO,
    0 otherwise. Can obtain negative examples by randomly sampling
    the genes that are not known to be associated with X

    # TODO: should we use the random split argument to the sklearn train/test split??

"""

# TODO: correctly compute AUC

import os
import numpy as np
from sklearn.cross_validation import train_test_split
from sklearn.cross_validation import StratifiedKFold
from sklearn import linear_model
from sklearn.metrics import roc_auc_score, roc_curve, auc
import argparse
import os
from math import log10
from random import shuffle

"""
def plot_roc(fprs, tprs, title):
    plt.figure()
    #for (fpr, tpr) in zip(fprs, tprs):
    plt.plot(fprs, tprs, 'gray')
    plt.plot([0, 1], [0, 1], 'k--')  # Plot the 50% line
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc="lower right")
    plt.show()
"""


def get_data(term, num_features=8555, cols=None):
    """
    Load and log-transform the data: x := log(x+1)

    :param term: The GO term
    :param num_features: # of features to extract (8555 is the full expression profile)
    :return: The data
    """

    if cols:
        num_features = len(cols)

    gene_features = np.empty((0, num_features))  # each row will be the feature profile for a given gene
    gene_ids = []

    # Get positive examples
    pos_file_name = '../data/experiment_inputs/' + term + '_pos.txt'
    pos_file = open(pos_file_name)
    for (i, line) in enumerate(pos_file):
        if i < 2:
            continue
        vals = line.rstrip().split('\t')
        gene_ids.append(vals[0])
        if cols:
            exp_levels = [vals[col+1] for col in cols]  # add 1 to skip over gene id column
        else:
            exp_levels = vals[1:]
        # Convert expression levels to log(level+1)
        exp_levels = [log10(float(exp_level)+1.0) for exp_level in exp_levels]
        gene_features = np.append(gene_features, [exp_levels], axis=0)
    pos_file.close()

    num_pos_examples = gene_features.shape[0]

    # Add on the negative examples
    neg_file_name = '../data/experiment_inputs/' + term + '_neg_0.txt'
    neg_file = open(neg_file_name)
    for (i, line) in enumerate(neg_file):
        if i < 2:
            continue
        vals = line.rstrip().split('\t')
        gene_ids.append(vals[0])
        if cols:
            exp_levels = [vals[col+1] for col in cols]  # add 1 to skip over gene id column
        else:
            exp_levels = vals[1:]
        # Convert expression levels to log(level+1)
        exp_levels = [log10(float(exp_level)+1.0) for exp_level in exp_levels]
        gene_features = np.append(gene_features, [exp_levels], axis=0)
    neg_file.close()

    num_neg_examples = gene_features.shape[0] - num_pos_examples

    labels = num_pos_examples * [1] + num_neg_examples * [0]

    # Normalize each feature (ie each column) to have 0 mean, unit variance
    for i in range(0, num_features):
        col = gene_features[:, i]
        col = col - np.mean(col)
        std_dev = np.std(col)
        if std_dev > 0:
            col = col / std_dev
        gene_features[:, i] = col

    indeces = range(0, gene_features.shape[0])
    # TODO: should I use random state?
    return gene_ids, train_test_split(gene_features, labels, indeces, test_size=0.33)  # , random_state=42)


def logistic_regresssion_L1(term, x_tr, x_te, y_tr, y_te, gene_ids_test, idx_te, server, tissues=None, rand_permute=False):
    num_folds = 3   # number of folds to use for cross-validation
    loss_function = 'l1'  # Loss function to use. Must be either 'l1' or 'l2'
    costs = np.logspace(-4, 4, 20)  # 10^(-start) to 10^stop in 10 logarithmic steps

    if rand_permute:
        # Randomly permute the labels
        shuffle(y_tr)
        shuffle(y_te)

    logreg_cv_L1 = linear_model.LogisticRegressionCV(Cs=costs, cv=num_folds, penalty=loss_function,
                                                     scoring='roc_auc', solver='liblinear', tol=0.0001)
    logreg_cv_L1.fit(x_tr, y_tr)
    best_c = logreg_cv_L1.C_
    pred_lr_cv_L1 = logreg_cv_L1.predict(x_te)
    conf_lr_cv_L1 = logreg_cv_L1.decision_function(x_te)
    prob_lr_cv_L1 = logreg_cv_L1.predict_proba(x_te)
    print_prediction_results('Cross-Validated Logistic Regression', logreg_cv_L1, y_te,
                             pred_lr_cv_L1, conf_lr_cv_L1, prob_lr_cv_L1,
                             other_info='Norm: ' + loss_function + ', # of Folds: ' + str(num_folds))
    if best_c != 0:
        print 'Best cost (after inverting to obtain true cost): ', 1.0 / best_c
    else:
        print 'Best cost (not inverted) is 0'

    # Save results
    if tissues:
        directory_path = 'results_tissue_specific/' + tissues[0] + '_' + str(server)
    else:
        directory_path = 'results_terms_with_tissues_v2_' + str(server)
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

    out_fname = directory_path + '/result_logreg_' + term + '.txt'
    save_prediction_results(out_fname, term, 'Logistic Regression with L1 Penalty', logreg_cv_L1,
                            y_te, pred_lr_cv_L1, conf_lr_cv_L1, prob_lr_cv_L1, gene_ids_test, num_folds,
                            tissue_set=tissues, costs=costs, best_cost=best_c)

    """
    copy_results_to_S3(out_fname, server)
    """


def copy_results_to_S3(out_fname, server):
    print 'Copying ', out_fname, ' back to S3.'
    os.system('aws s3 cp ' + out_fname + ' s3://stanfordgtex/GO_Prediction/GO_Prediction_Results/' + out_fname)


def copy_zipped_results_to_S3(server):
    # Create tarball
    os.system('tar -zcvf results_server_' + str(server) + '.tar.gz results_server_' + str(server))

    # Send the whole compressed directory to S3
    os.system('aws s3 cp results_server_' + str(server) + '.tar.gz ' +
              's3://stanfordgtex/GO_Prediction/GO_Prediction_Results/results_server_' + str(server) + '.tar.gz')


def print_prediction_results(model, fit, labels, predictions, conf_scores, pred_prob, other_info=None):
    print 20*'-'
    print model
    print 20*'-'
    prob_of_pos = [p[1] for p in pred_prob]
    print 'ROC AUC Score: ', roc_auc_score(labels, prob_of_pos)
    print 'Pred prob: ', pred_prob
    print 'Conf scores: ', conf_scores
    print labels
    print predictions
    print 'Number of nonzero coefficients: ', np.count_nonzero(fit.coef_)
    print 'Total number of coefficients: ', len(fit.coef_[0])

    if other_info:
        print other_info

    return


def save_prediction_results(fname, GO_id, model, fit, test_labels, preds, conf_scores, prob_preds, gene_ids_test, n_folds=None,
                            tissue_set=None, costs=None, best_cost=None, other_info=None):
    """
    Save results of prediction for a given GO term to a text file.

    :param fname: Name of output file
    :param GO_id: ID of the GO term
    :param model: Type of model used (e.g. "Logistic Regression with L1 Norm")
    :param test_labels Labels of the test data
    :param preds: Binary predictions for test data
    :param gene_ids_test: Gene ids of the genes in the test set
    :param n_folds: # of folds used for cross-validation
    :param tissue_set: Set of tissues used. If no option is specified, it's assumed that all tissues were used
    :param costs: If applicable, range of costs used
    :param best_cost: If applicable, the best cost parameter (as determined by cross-validation)
    :param other_info: Other info to store in header of output file
    :return: None
    """

    out_file = open(fname, 'w')
    out_file.write('# Prediction results for GO term: ' + GO_id + '\n')
    out_file.write('# Model used: ' + model + '\n')
    prob_of_pos = [p[1] for p in prob_preds]
    auc_score = roc_auc_score(test_labels, prob_of_pos)

    out_file.write('# ROC AUC Score: ' + str(auc_score) + '\n')
    if tissue_set:
        out_file.write('# Tissues used: ' + str(tissue_set) + '\n')
    else:
        out_file.write('# All tissues were included\n')
    if n_folds:
        out_file.write('# Number of folds used for cross-validation: ' + str(n_folds) + '\n')
    if costs is not None:
        out_file.write('# Range of cost parameters:\n')
        for c in costs[:-1]:
            out_file.write(str(c) + '\t')
        out_file.write(str(costs[-1]) + '\n')
    if best_cost:
        out_file.write('# Best cost parameter (determined by CV): ' + str(best_cost) + '\n')
    if other_info:
        out_file.write('# ' + other_info + '\n')

    # Write out the coefficients
    out_file.write('# Coefficients:\n')
    coefficients = fit.coef_[0]
    for c in coefficients[:-1]:
        out_file.write(str(c) + '\t')
    out_file.write(str(coefficients[-1]) + '\n')

    out_file.write('# Gene_ID\tLabel\tPrediction\tDecision_Func_Score\tProb_of_Pos\n')
    for id,label,pred,conf,prob in zip(gene_ids_test, test_labels, preds, conf_scores, prob_preds):
        out_file.write(id + '\t' + str(label) + '\t' + str(pred) + '\t' + str(conf) + '\t' + str(prob[1]) + '\n')

    out_file.close()


def get_go_terms():
    f_name = '../data/GO_terms_final_gene_counts.txt'
    GO_counts_file = open(f_name)

    terms = []
    for (i, line) in enumerate(GO_counts_file):
        if i < 2:
            continue
        term = line.split('\t')[0]
        terms.append(term)

    return terms


def get_tissue_list(tissue_fpath):
    tissue_file = open(tissue_fpath)
    for line in tissue_file:
        tissues = line.rstrip().split('\t')
        break
    return tissues


def get_tissues_to_cols(tissue_list):
    tissues_to_cols = {}
    for tissue in tissue_list:
        cols = []
        meta_fname = '../data/tissue_metadata/tissue_meta_' + tissue + '.txt'
        meta_file = open(meta_fname)
        for (i, line) in enumerate(meta_file):
            if i < 1:
                continue
            col = int(line.split('\t')[0])
            cols.append(col)
        tissues_to_cols[tissue] = cols
    return tissues_to_cols


def run_prediction(GO_term, server_no, tissues=None, cols=None):
    genes_train_test, train_test = get_data(GO_term, cols=cols)
    X_train, X_test, y_train, y_test, idx_train, idx_test = train_test
    # Arrange gene ids to match ordering of test set
    test_genes = [genes_train_test[idx] for idx in idx_test]

    # Run Logistic Regression
    logistic_regresssion_L1(GO_term, X_train, X_test, y_train, y_test, test_genes, idx_test, server_no, tissues=tissues)


def predict_GO_to_genes(n_servers, server_no, tissue_specific=False, tissue_fpath=None):

    if tissue_specific:
        # Determine which columns each tissue corresponds to
        tissues = get_tissue_list(tissue_fpath)
        tissues_to_cols = get_tissues_to_cols(tissues)

    GO_terms = get_go_terms()
    GO_terms.reverse()   # Process GO terms in order from least # of genes to most

    terms_to_process = len(GO_terms)/n_servers
    print 'This program will process approximately ', terms_to_process, ' GO terms'

    # todo: 670-720 for the unscaled results (both rand and non rand for the log transformed data)
    for (idx, GO_term) in enumerate(GO_terms):

        if idx % n_servers != server_no:
            # Not for this server
            continue

        print '=' * 30
        print '=' * 30
        print idx, 'th GO term: ', GO_term
        print 'Approximately ', (terms_to_process-idx/n_servers), ' terms left to process'

        if tissue_specific:
            for tissue in tissues:
                run_prediction(GO_term, server_no, tissues=[tissue], cols=tissues_to_cols[tissue])
                break
        else:
            run_prediction(GO_term, server_no)
            """
            genes_train_test, train_test = get_data(GO_term, tissues)
            X_train, X_test, y_train, y_test, idx_train, idx_test = train_test
            # Arrange gene ids to match ordering of test set
            test_genes = [genes_train_test[idx] for idx in idx_test]

            # Run Logistic Regression
            logistic_regresssion_L1(GO_term, X_train, X_test, y_train, y_test, test_genes, idx_test, server_no)
            """
        break

"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--n_servers', default=1, type=int, required=False,
                        help='Tells the program how many servers you\'re using')
    parser.add_argument('--server_no', default=0, type=int, required=False,
                        help='Which server # this program is running on')
    args = parser.parse_args()
    num_servers = args.n_servers
    server_number = args.server_no
    print '# of servers to be used: ', str(num_servers)
    print 'This program is running on server: ', str(server_number)

    # Make predictions using expression across all tissues
    #predict_GO_to_genes(num_servers, server_number)

    # Make predictions using expression in individual tissues
    predict_GO_to_genes(num_servers, server_number, tissue_specific=True, tissue_fpath='../data/tissues.txt')





