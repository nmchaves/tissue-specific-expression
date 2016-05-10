"""
    Consider biological process X
    Training example: gene
    Label (binary): 1 if gene is associated with X in GO,
    0 otherwise. Can obtain negative examples by randomly sampling
    the genes that are not known to be associated with X

    # TODO: cross-validation, AUC score, clean up pipeline, add in tissue
    selection to the pipeline

    # TODO: use only 53 tissues
"""

import numpy as np
from sklearn import linear_model
from sklearn.svm import SVC
import GO_utils
import utils


"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    gene2go_file_path = '../data/gene2go.txt' # If file doesn't exist, then run gene2go = download_ncbi_associations()
    rpkm_file_path = '../../CS341_Data/transcript_rpkm_in_go_nonzero_exp.txt'
    gene_count_file_path = '../data/supp_GO_term_gene_counts.txt'
    biomart_file_path = '../data/biomart_ensembl_to_entrez.txt'
    #gene2go_file_path = '../local_data/gene2go.txt' # If file doesn't exist, then run gene2go = download_ncbi_associations()
    #rpkm_file_path = '../local_data/transcript_rpkm_in_go.txt'
    sample_tissue_path = '../data/sampleID_tissue.txt'

    GO_PROCESS_IDs = GO_utils.get_go_terms(biomart_file_path, gene2go_file_path, gene_count_file_path, top=10)

    # GO:0007596
    (GO_PROCESS_ID, ensembl_ids) = GO_PROCESS_IDs[7]
    print GO_PROCESS_ID

    go_evidence_codes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']
    #ensembl_ids = get_ensembl_ids(GO_PROCESS_ID, biomart_file_path, ev_codes=go_evidence_codes)

    # NUM_FEATURES = 8555
    NUM_FEATURES = 53

    # 1st Pass Through Dataset: Obtain positive training examples
    gene_features, positive_example_rows, gene_ids_ordered, num_transcripts = \
        GO_utils.get_positive_examples(rpkm_file_path, sample_tissue_path, ensembl_ids)

    print 'After pass 1 (inserting positive examples), gene feature matrix has dimension: ', gene_features.shape
    num_positive_examples = len(positive_example_rows)
    num_negative_examples = num_positive_examples
    num_examples = num_positive_examples + num_negative_examples

    # 2nd Pass through dataset: Obtain an equal number of negative training exmaples
    neg_rows = utils.rand_sample_exclude(range(0, num_transcripts), num_positive_examples, exclude=positive_example_rows)

    gene_features_neg, gene_ids_ordered_neg = \
        GO_utils.get_negative_examples(rpkm_file_path, sample_tissue_path, neg_rows)
    gene_features = np.append(gene_features, gene_features_neg, axis=0)
    gene_ids_ordered += gene_ids_ordered_neg

    print 'After pass 2 (inserting negative examples), gene feature matrix has dimension: ', gene_features.shape

    # Vector of labels for each example
    labels = num_positive_examples * [1] + num_negative_examples * [0]

    train_data, test_data = utils.split_data(gene_features, labels, gene_ids_ordered, train_set_size=0.7)

    # TODO: train_data.gene_features
    '''
    logreg = linear_model.LogisticRegression(C=1e5)
    logreg.fit(gene_features_train, labels_train)
    pred_lr = logreg.predict(gene_features_test)
    print_prediction_results('Logistic Regression', labels_test, pred_lr)
    '''

    # Logistic Regression with Cross-Validation, L1 Norm (must use liblinear solver for L1)
    #costs = []
    '''
    num_folds = 5   # number of folds to use for cross-validation
    loss_function = 'l1'  # Loss function to use. Must be either 'l1' or 'l2'
    logreg_cv_L1 = linear_model.LogisticRegressionCV(cv=num_folds, penalty=loss_function, solver='liblinear')
    logreg_cv_L1.fit(gene_features_train, labels_train)
    pred_lr_cv_L1 = logreg_cv_L1.predict(gene_features_test)
    print_prediction_results('Cross-Validated Logistic Regression', labels_test, pred_lr_cv_L1,
                             other_info='Norm: ' + loss_function + ', # of Folds: ' + str(num_folds))
    '''

    num_folds = 5   # number of folds to use for cross-validation
    loss_function = 'l2'  # Loss function to use. Must be either 'l1' or 'l2'
    logreg_cv_L2 = linear_model.LogisticRegressionCV(cv=num_folds, penalty=loss_function)
    logreg_cv_L2.fit(gene_features_train, labels_train)
    pred_lr_cv_L2 = logreg_cv_L2.predict(gene_features_test)
    utils.print_prediction_results('Cross-Validated Logistic Regression', labels_test, pred_lr_cv_L2
                             , other_info='Norm: ' + loss_function + ', # of Folds: ' + str(num_folds))

    # SVM
    clf = SVC()
    clf.fit(gene_features_train, labels_train)

    pred_svm = clf.predict(gene_features_test)
    utils.print_prediction_results('SVM', labels_test, pred_svm)











