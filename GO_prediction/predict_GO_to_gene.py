#! /home/ec2-user/anaconda2/bin/python
"""
    Script for running prediction tasks. In particular, this script
    can be used to predict whether a given GO term is associated with
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
from sklearn.cross_validation import StratifiedKFold
from sklearn import linear_model
from sklearn.metrics import roc_auc_score, roc_curve, auc
import argparse
import os
from math import log10
from random import shuffle


def get_tissue_list(tissue_fpath):
    tissue_file = open(tissue_fpath)
    for line in tissue_file:
        tissues = line.rstrip().split('\t')
        break
    return tissues


def get_data(term, num_features=8555, cols=None, pca_dimen=None, median=False):
    """
    Load and log-transform the data: x := log(x+1)

    :param term: The GO term
    :param num_features: # of features to extract (8555 is the full expression profile)
    :return: The data
    """

    if cols:
        num_features = len(cols)
    elif pca_dimen:
        num_features = pca_dimen * NUM_TISSUES
    elif median:
        num_features = NUM_TISSUES

    gene_features = np.empty((0, num_features))  # each row will be the feature profile for a given gene
    gene_ids = []

    # Get positive examples
    pos_file_name = exp_input_file_path + term + '_pos.txt'
    pos_file = open(pos_file_name)
    for (i, line) in enumerate(pos_file):
        if i < 2:
            continue
        vals = line.rstrip().split('\t')[1:]  # Skip over the 1st column, which is a row index
        gene_ids.append(vals[0])
        if cols:
            exp_levels = [float(vals[col+1]) for col in cols]  # add 1 to skip over gene id column
        else:
            exp_levels = [float(val) for val in vals[1:]]

        if not pca_dimen:
            # Convert expression levels to log(level+1)
            exp_levels = [log10(float(exp_level)+1.0) for exp_level in exp_levels]
        gene_features = np.append(gene_features, [exp_levels], axis=0)
    pos_file.close()

    num_pos_examples = gene_features.shape[0]

    # Add on the negative examples
    neg_file_name = exp_input_file_path + term + '_neg_' + str(neg_set) + '.txt'
    neg_file = open(neg_file_name)
    for (i, line) in enumerate(neg_file):
        if i < 2:
            continue
        vals = line.rstrip().split('\t')[1:]
        gene_ids.append(vals[0])
        if cols:
            exp_levels = [float(vals[col+1]) for col in cols]  # add 1 to skip over gene id column
        else:
            exp_levels = [float(val) for val in vals[1:]]

        if not pca_dimen:
            # Convert expression levels to log(level+1)
            exp_levels = [log10(float(exp_level)+1.0) for exp_level in exp_levels]
        gene_features = np.append(gene_features, [exp_levels], axis=0)
    neg_file.close()

    num_neg_examples = gene_features.shape[0] - num_pos_examples

    labels = num_pos_examples * [1] + num_neg_examples * [0]

    if not pca_dimen:
        # Normalize each feature (ie each column) to have 0 mean, unit variance
        for i in range(0, num_features):
            col = gene_features[:, i]
            col = col - np.mean(col)
            std_dev = np.std(col)
            if std_dev > 0:
                col = col / std_dev
            gene_features[:, i] = col

    indeces = range(0, gene_features.shape[0])
    # TODO: should I use random state? (good for reproducible results when debugging)
    return gene_ids, train_test_split(gene_features, labels, indeces, stratify=labels, test_size=0.33)  # , random_state=42)


def logistic_regresssion(term, x_tr, x_te, y_tr, y_te, gene_ids_test, idx_te, server,
                         loss_function='l2', tissues=None, rand_permute=False, pca_dimen=None):
    num_folds = 5   # number of folds to use for cross-validation
    #loss_function = 'l1'  # Loss function to use. Must be either 'l1' or 'l2'
    costs = np.logspace(-4, 4, 20)  # 10^(-start) to 10^stop in 10 logarithmic steps

    if rand_permute:
        # Randomly permute the labels
        shuffle(y_tr)
        shuffle(y_te)

    logreg_cv = linear_model.LogisticRegressionCV(Cs=costs, cv=num_folds, penalty=loss_function,
                                                  scoring='roc_auc', solver='liblinear', tol=0.0001)
    logreg_cv.fit(x_tr, y_tr)
    best_c = logreg_cv.C_
    pred_lr_cv = logreg_cv.predict(x_te)
    conf_lr_cv = logreg_cv.decision_function(x_te)
    prob_lr_cv = logreg_cv.predict_proba(x_te)

    if not AWS:
        print_prediction_results('Cross-Validated Logistic Regression', logreg_cv, y_te,
                                 pred_lr_cv, conf_lr_cv, prob_lr_cv,
                                 other_info='Norm: ' + loss_function + ', # of Folds: ' + str(num_folds))
        if best_c != 0:
            print 'Best cost (after inverting to obtain true cost): ', 1.0 / best_c
        else:
            print 'Best cost (not inverted) is 0'

    # Save results
    directory_path = output_dir
    if tissues:
        directory_path += '/' + term + '_neg_' + str(neg_set)
        out_fname = directory_path + '/logreg_' + tissues[0] + '.txt'
    else:
        out_fname = directory_path + '/logreg_' + term + '.txt'
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

    if loss_function == 'l1':
        model = 'Logistic Regression with L1 Penalty'
    elif loss_function == 'l2':
        model = 'Logistic Regression with L2 Penalty'
    else:
        model = 'Logistic Regression. Penalty is neither l1 nor l2.'

    save_prediction_results(out_fname, term, model, logreg_cv,
                            y_te, pred_lr_cv, conf_lr_cv, prob_lr_cv, gene_ids_test, num_folds,
                            tissue_set=tissues, costs=costs, best_cost=best_c, pca_dimen=pca_dimen)

    if AWS:
        copy_results_to_S3(out_fname, server)



def copy_results_to_S3(out_fname, server):

    # Before saving, convert characters to their escaped versions
    out_fname = repr(out_fname)

    print 'Copying ', out_fname, ' back to S3.'
    os.system('aws s3 cp ' + out_fname + ' s3://stanfordgtex/GO_Prediction/Median_GO_Prediction_Results/' + out_fname)


def copy_zipped_results_to_S3(server, tissue_specific=False):
    # Create tarball
    tar_name = 'full_results_'
    if tissue_specific:
        tar_name += '1_tissue_'
    else:
        tar_name += 'all_tissues_'
    tar_name += 'loss_' + loss + '_neg_' + str(neg_set) + '_server_' + str(server)
    os.system('tar -zcvf ' + tar_name + '.tar.gz ' + tar_name)

    # Send the whole compressed directory to S3
    os.system('aws s3 cp ' + tar_name + '.tar.gz' +
              ' s3://stanfordgtex/GO_Prediction/Median_GO_Prediction_Results/' + tar_name + '.tar.gz')


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
                            tissue_set=None, costs=None, best_cost=None, other_info=None, pca_dimen=None):
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
        out_file.write('# Tissues used:\n')
        for tiss in tissue_set[:-1]:
            out_file.write(str(tiss) + '\t')
        out_file.write(str(tissue_set[-1] + '\n'))
    else:
        out_file.write('# All tissues were included\n')
    if pca_dimen:
        out_file.write('# PCA used on features. Kept top ' + str(pca_dimen) + ' components.\n')
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
    if AWS:
        f_name = 'GO_terms_final_gene_counts.txt'
    else:
        f_name = '../data/GO_terms_final_gene_counts.txt'
    GO_counts_file = open(f_name)

    terms = []
    for (i, line) in enumerate(GO_counts_file):
        if i < 2:
            continue
        term = line.split('\t')[0]
        terms.append(term)

    return terms


def get_tissues_to_cols(tissue_list, pca_dimen=None, median=False):
    tissues_to_cols = {}

    if pca_dimen:
        for (i, tissue) in enumerate(tissue_list):
            col_start = pca_dimen * i
            col_stop = col_start + pca_dimen
            tissues_to_cols[tissue] = range(col_start, col_stop)
    elif median:
        for (i, tissue) in enumerate(tissue_list):
            tissues_to_cols[tissue] = [i]
    else:
        for tissue in tissue_list:
            cols = []
            if AWS:
                meta_fname = 'tissue_metadata/tissue_meta_' + tissue + '.txt'
            else:
                meta_fname = '../data/tissue_metadata/tissue_meta_' + tissue + '.txt'
            meta_file = open(meta_fname)
            for (i, line) in enumerate(meta_file):
                if i < 1:
                    continue
                col = int(line.split('\t')[0])
                cols.append(col)
            tissues_to_cols[tissue] = cols
    return tissues_to_cols


def run_prediction(GO_term, server_no, tissues=None, cols=None, loss=None, pca_dimen=None, median=None):
    genes_train_test, train_test = get_data(GO_term, cols=cols, pca_dimen=pca_dimen, median=median)
    X_train, X_test, y_train, y_test, idx_train, idx_test = train_test

    # Arrange gene ids to match ordering of test set
    test_genes = [genes_train_test[idx] for idx in idx_test]

    # Run Logistic Regression
    logistic_regresssion(GO_term, X_train, X_test, y_train, y_test, test_genes, idx_test,
                         server_no, loss_function=loss, tissues=tissues)


def predict_GO_to_genes(n_servers, server_no, tissue_specific=False, tissue_fpath=None,
                        loss=None, pca_dimen=None, median=False):

    if tissue_specific:
        # Determine which columns each tissue corresponds to
        tissues = get_tissue_list(tissue_fpath)
        tissues_to_cols = get_tissues_to_cols(tissues, pca_dimen=pca_dimen, median=median)

    GO_terms = get_go_terms()
    GO_terms.reverse()   # Process GO terms in order from least # of genes to most

    terms_to_process = len(GO_terms)/n_servers
    print 'This program will process approximately ', terms_to_process, ' GO terms'

    for (idx, GO_term) in enumerate(GO_terms):

        if idx % n_servers != server_no:
            # Not for this server
            continue

        print '=' * 30
        print '=' * 30
        print '=' * 30
        print '=' * 30
        print idx, 'th GO term: ', GO_term
        print 'Approximately ', (terms_to_process-idx/n_servers), ' terms left to process'

        if tissue_specific:
            for tissue in tissues:
                run_prediction(GO_term, server_no, tissues=[tissue], cols=tissues_to_cols[tissue],
                               loss=loss, pca_dimen=pca_dimen, median=median)
        else:
            run_prediction(GO_term, server_no, loss=loss, pca_dimen=pca_dimen, median=median)

    if AWS:
        copy_zipped_results_to_S3(server_no, tissue_specific)

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
    parser.add_argument('--neg_set', default=0, type=int, required=False,
                        help='Which negative set to used.')
    parser.add_argument('--aws', default=True, type=bool, required=False,
                        help='Whether or not the code is being run on AWS.')
    parser.add_argument('--single_tissue', default=True, type=bool, required=False,
                        help='If True, each prediction problem only used features from 1 tissue')

    args = parser.parse_args()
    num_servers = args.n_servers
    server_number = args.server_no
    neg_set = args.neg_set
    AWS = True #args.aws
    single_tissue = True #args.single_tissue

    loss = 'l2'
    pca_dimen = None
    median = True

    if AWS:
        exp_input_file_path = 'median_experiment_inputs/'
        tissue_path = 'tissues.txt'
    else:
        exp_input_file_path = '../../CS341_Data/median_experiment_inputs/'
        tissue_path = '../data/tissues.txt'

    NUM_TISSUES = 53

    print '# of servers to be used: ', str(num_servers)
    print 'This program is running on server: ', str(server_number)
    print 'Obtaining experiment inputs from negative set: ', str(neg_set)

    if pca_dimen:
        output_dir = 'pca_results_'
    elif median:
        output_dir = 'median_results_'
    else:
        output_dir = 'full_results_'

    if single_tissue:
        output_dir += '1_tissue_loss_' + loss + '_neg_' + str(neg_set) + '_server_' + str(server_number)

        # Make predictions using expression in individual tissues
        predict_GO_to_genes(num_servers, server_number, tissue_specific=True, tissue_fpath=tissue_path, loss=loss,
                            pca_dimen=pca_dimen, median=median)
    else:
        output_dir += 'all_tissues_loss_' + loss + '_neg_' + str(neg_set) + '_server_' + str(server_number)

        # Make predictions using expression across all tissues
        predict_GO_to_genes(num_servers, server_number, loss=loss, pca_dimen=pca_dimen, median=median)

