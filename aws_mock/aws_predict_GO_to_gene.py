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

import numpy as np
from sklearn.cross_validation import train_test_split
from sklearn import linear_model
from sklearn.metrics import roc_auc_score
import argparse
import os


def get_data(term, num_features=8555):

    gene_features = np.empty((0, num_features))  # each row will be the feature profile for a given gene
    gene_ids = []

    # Get positive examples
    pos_file_name = 'experiment_inputs_subset/' + term + '_pos.txt'
    pos_file = open(pos_file_name)
    skipLines = 2
    for (i, line) in enumerate(pos_file):
        if i < skipLines:
            continue
        vals = line.rstrip().split('\t')
        gene_ids.append(vals[0])
        exp_levels = vals[1:]
        exp_levels = [float(exp_level) for exp_level in exp_levels]
        gene_features = np.append(gene_features, [exp_levels], axis=0)
    pos_file.close()

    num_pos_examples = gene_features.shape[0]

    # Add on the negative examples
    neg_file_name = 'experiment_inputs_subset/' + term + '_neg_0.txt'
    neg_file = open(neg_file_name)
    skipLines = 2
    for (i, line) in enumerate(neg_file):
        if i < skipLines:
            continue
        vals = line.rstrip().split('\t')
        gene_ids.append(vals[0])
        exp_levels = vals[1:]
        exp_levels = [float(exp_level) for exp_level in exp_levels]
        gene_features = np.append(gene_features, [exp_levels], axis=0)
    pos_file.close()

    num_neg_examples = gene_features.shape[0] - num_pos_examples

    labels = num_pos_examples * [1] + num_neg_examples * [0]

    indeces = range(0, gene_features.shape[0])
    # TODO: should I use random state?
    return gene_ids, train_test_split(gene_features, labels, indeces, test_size=0.33)  # , random_state=42)


def logistic_regresssion_L1(term, x_tr, x_te, y_tr, y_te, gene_ids_test, idx_te, server):
    num_folds = 3   # number of folds to use for cross-validation
    loss_function = 'l1'  # Loss function to use. Must be either 'l1' or 'l2'
    costs = np.logspace(-4, 4, 20)  # 10^(-start) to 10^stop in 10 logarithmic steps
    logreg_cv_L1 = linear_model.LogisticRegressionCV(Cs=costs, cv=num_folds, penalty=loss_function,
                                                     scoring='roc_auc', solver='liblinear', tol=0.0001)
    logreg_cv_L1.fit(x_tr, y_tr)
    best_c = logreg_cv_L1.C_
    pred_lr_cv_L1 = logreg_cv_L1.predict(x_te)
    print_prediction_results('Cross-Validated Logistic Regression', y_te, pred_lr_cv_L1,
                             other_info='Norm: ' + loss_function + ', # of Folds: ' + str(num_folds))
    if best_c != 0:
        print 'Best cost (after inverting to obtain true cost): ', 1.0 / best_c
    else:
        print 'Best cost (not inverted) is 0'

    # Save results
    directory_path = 'results_server_' + str(server)
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
    out_fname = directory_path + '/new_result_logreg_' + term + '.txt'
    save_prediction_results(out_fname, term, 'Logistic Regresion with L1 Penalty',
                            y_te, pred_lr_cv_L1, gene_ids_test, num_folds,
                            tissue_set=None, costs=costs, best_cost=best_c)


def print_prediction_results(model, labels, predictions, other_info=None):
    print 20*'-'
    print model
    print 20*'-'
    print 'ROC AUC Score: ', roc_auc_score(labels, predictions)

    if other_info:
        print other_info
    return


def save_prediction_results(fname, GO_id, model, test_labels, preds, gene_ids_test, n_folds=None,
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
    auc_score = roc_auc_score(test_labels, preds)

    out_file.write('# ROC AUC Score: ' + str(auc_score) + '\n')
    if tissue_set:
        out_file.write('# Tissues used: ' + tissue_set + '\n')
    else:
        out_file.write('# All tissues were included\n')
    if n_folds:
        out_file.write('# Number of folds used for cross-validation: ' + str(n_folds) + '\n')
    if costs is not None:
        out_file.write('# Range of cost parameters: ' + str(costs) + '\n')
    if best_cost:
        out_file.write('# Best cost parameter (determined by CV): ' + str(best_cost) + '\n')
    if other_info:
        out_file.write('# ' + other_info + '\n')

    out_file.write('# Gene ID\tLabel\tPrediction\n')
    for id,label,pred in zip(gene_ids_test, test_labels, preds):
        out_file.write(id + '\t' + str(label) + '\t' + str(pred) + '\n')

    out_file.close()


def get_go_terms():
    f_name = 'GO_terms_final_gene_counts_subset.txt'
    GO_counts_file = open(f_name)

    terms = []
    for (i, line) in enumerate(GO_counts_file):
        if i < 2:
            continue
        term = line.split('\t')[0]
        terms.append(term)

    return terms


"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--n_servers', default=3, type=int, required=False,
                    help='Tells the program how many servers you\'re using')
    parser.add_argument('--server_no', default=0, type=int, required=False,
                        help='Which server # this program is running on')
    args = parser.parse_args()
    num_servers = args.n_servers
    server_no = args.server_no
    print '# of servers to be used: ', str(num_servers)
    print 'This program is running on server: ', str(server_no)

    # TODO: exceptions for invalid server #'s, etc.

    GO_terms = get_go_terms()
    GO_terms.reverse()   # Process GO terms in order from least # of genes to most

    print 'This program will process approximately ', len(GO_terms)/num_servers, ' GO terms'

    for (idx, GO_term) in enumerate(GO_terms):

        if idx % num_servers != server_no:
            # Not for this server
            continue

        print GO_term

        genes_train_test, train_test = get_data(GO_term)
        X_train, X_test, y_train, y_test, idx_train, idx_test = train_test
        # Arrange gene ids to match ordering of test set
        test_genes = [genes_train_test[idx] for idx in idx_test]

        print idx_test
        print y_train
        print y_test
        print test_genes

        # Run Logistic Regression
        logistic_regresssion_L1(GO_term, X_train, X_test, y_train, y_test, test_genes, idx_test, server_no)
