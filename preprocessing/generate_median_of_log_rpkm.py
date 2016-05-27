"""
    This file goes through all GO terms of interest and creates the positive and negative
    examples for running a prediction problem. A positive example is a gene that is known
    to be associated with the GO term, while a negative example is a randomly sampled gene
    from the set of genes that are not known to be associated with the GO term.

"""

from GO_Evidence_Codes import EvidenceCodes
from os import remove
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

import sys
# sys.path.insert(0, '../GO_prediction')
# import GO_utils
# import utils
# from utils import get_tissue_list

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

def normalize_columns(mtx_in):
    mtx_out = np.zeros(mtx_in.shape)
    num_features = mtx_in.shape[1]
    for i in range(num_features):
        col = mtx_in[:, i]
        col = col - np.mean(col)
        std_dev = np.std(col)
        if std_dev > 0:
            col = col / std_dev
        mtx_out[:, i] = col
    return mtx_out

'''
*********************
        Main
*********************
'''

if __name__ == "__main__":

    # --------------------------------------
    # input and output files
    # --------------------------------------
    # rpkm_file_path = '../data/local_large/transcript_rpkm_in_go_nonzero_exp.txt'
    # rpkm_file_out = '../data/local_large/log_norm_pca_transcript_rpkm_in_go_nonzero_exp.txt'
    rpkm_file_path = '../data/small_example_data/small_transcript_rpkm_in_go_nonzero_exp.txt'
    rpkm_file_out = '../data/small_example_data/small_log_median_transcript_rpkm_in_go_nonzero_exp.txt'

    # --------------------------------------
    # load tissue information
    # --------------------------------------
    tissues = get_tissue_list('../data/tissues.txt');
    tissues_to_cols = get_tissues_to_cols(tissues)
    num_features = 0
    for tissue in tissues_to_cols:
        num_features += len(tissues_to_cols[tissue])
    print 'Loaded '+str(num_features)+' features. Number of tissues:' + str(len(tissues))

    # --------------------------------------
    # load count data and meta information
    # --------------------------------------
    n_header_lines = 1
    n_meta_fields = 4
    # 1. count number of genes
    rpkm_file = open(rpkm_file_path)
    num_genes = 0
    for (i, line) in enumerate(rpkm_file):
        if i >= n_header_lines:
            num_genes += 1
    rpkm_file.close()
    print 'Reading: ' + rpkm_file_path
    print 'Number of genes: ' + str(num_genes)
    # 2. read the full rpkm matrix
    full_log_mtx = np.empty((num_genes, num_features))  # each row will be the feature profile for a given gene
    gene_info = ["NULL"]*num_genes
    print 'Reading:' + rpkm_file_path
    rpkm_file = open(rpkm_file_path)
    for (i, line) in enumerate(rpkm_file):
        vals = line.rstrip().split('\t')
        if i < n_header_lines: 
            first_fields = '\t'.join(vals[0:(n_meta_fields-1)]) # print header later
            continue
        idx = i - n_header_lines  # gene index 
        if (idx % 1000 == 0):
            print '    ' + str(idx) + ' genes read'
        gene_info[idx] = '\t'.join(vals[0:(n_meta_fields-1)])
        exp_levels = vals[n_meta_fields:]
        exp_levels = [np.log10(float(exp_level)+1.0) for exp_level in exp_levels] # log transform
        full_log_mtx[idx,:] = exp_levels 
    rpkm_file.close()
    print 'Completed reading full matrix. Dimensions: ' + str(full_log_mtx.shape)

    # --------------------------------------
    # calculate median for each tissue
    # --------------------------------------
    reduced_mtx = np.zeros([full_log_mtx.shape[0],len(tissues)]);
    reduced_cols= tissues 
    for idx, tissue in enumerate(tissues):
        print 'Taking median of ' + tissue 
        tissue_mtx = full_log_mtx[:,tissues_to_cols[tissue]] 
        reduced_mtx[:,idx] = np.median(tissue_mtx,axis=1)
    print 'Reduced matrix size: ', str(reduced_mtx.shape)

    # --------------------------------------
    # save reduced matrix to file 
    # --------------------------------------
    print 'Saving to: ' + rpkm_file_out
    with open(rpkm_file_out , mode='wt') as myfile: 
        myfile.write(first_fields)
        myfile.write('\t')
        myfile.write('\t'.join(reduced_cols))
        myfile.write('\n')
        for idx in range(len(gene_info)):
            myfile.write(gene_info[idx])
            myfile.write('\t')
            myfile.write('\t'.join(map(str,reduced_mtx[idx,:])))
            myfile.write('\n')

