"""
    File for predicting whether a given GO term is associated with
    various genes.

    Consider biological process X
    Training example: gene
    Label (binary): 1 if gene is associated with X in GO,
    0 otherwise. Can obtain negative examples by randomly sampling
    the genes that are not known to be associated with X

    # TODO: It would be best to use the random split argument to the sklearn train/test split??

"""

import os
import numpy as np
from sklearn.cross_validation import train_test_split
from sklearn.cross_validation import StratifiedKFold
from sklearn.decomposition import PCA
from sklearn import linear_model
from sklearn.metrics import roc_auc_score, roc_curve, auc
import argparse
import os
from math import log10
from random import shuffle
from utils import get_tissue_list

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


def get_mtx(term, fpath, num_features=8555, cols=None):
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
    # pos_file_name = '../../CS341_Data/experiment_inputs/' + term + '_pos.txt'
    pos_file_name = fpath + term + '_pos.txt'
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
    # neg_file_name = '../../CS341_Data/experiment_inputs/' + term + '_neg_0.txt'
    neg_file_name = fpath + term + '_neg_0.txt'
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
    return gene_features, gene_ids, labels 
    # return gene_ids, train_test_split(gene_features, labels, indeces, stratify=labels, test_size=0.33)  # , random_state=42)


def get_reduced_data(term, fpath, n_pcomp = 5):
    print '*Taking the top '+ str(n_pcomp) + ' principle components of each tissue'
    # extract tissue meta information
    tissues = get_tissue_list('../data/tissues.txt');
    tissues_to_cols = get_tissues_to_cols(tissues)
    # load an example matrix
    full_mtx, gene_ids, labels = get_mtx(term,fpath, num_features=8555)
    
    # perform pca for each tissue
    pca_tissues_to_cols = {}
    reduced_mtx = np.zeros([full_mtx.shape[0],n_pcomp*len(tissues)]);
    for idx, tissue in enumerate(tissues):
        pca_tissues_to_cols[tissue] = range(n_pcomp*idx, n_pcomp*(idx+1))
        tissue_mtx = full_mtx[:,tissues_to_cols[tissue]]
        assert(tissue_mtx.shape[1] >= n_pcomp)
        pca = PCA(n_components=n_pcomp)
        reduced_tissue_mtx = pca.fit_transform(tissue_mtx)
        reduced_mtx[:,pca_tissues_to_cols[tissue]] = reduced_tissue_mtx
        # print tissue_mtx.shape
        # print reduced_tissue_mtx.shape
        # print pca_tissues_to_cols[tissue]
    print '*Data matrix size: ', str(reduced_mtx.shape)
    # create train test sets from the dimension reduced matrix
    indeces = range(0, reduced_mtx.shape[0])
    train_test = train_test_split(reduced_mtx, labels, indeces, stratify=labels, test_size=0.33)  # , random_state=42)
    return gene_ids, train_test


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

"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    gene_ids,train_test = get_reduced_data('GO-0000578','../data/small_example_data/', n_pcomp=5)

    # print gene_ids
    # print train_test 
    

    
