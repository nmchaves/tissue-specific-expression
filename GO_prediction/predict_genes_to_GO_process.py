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
from sklearn.grid_search import GridSearchCV

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
    obo_file_path = '../data/go-basic.obo'

    # GO Evidence Codes
    exp_ev_codes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']
    comp_an_ev_codes = ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD', 'RCA']
    auth_stmt_ev_codes = ['TAS', 'NAS']
    cur_ev_codes = ['IC', 'ND']
    elec_ev_codes = ['IEA']


    NUM_FEATURES = 8555
    #NUM_FEATURES = 53

    ev_codes = exp_ev_codes + comp_an_ev_codes + auth_stmt_ev_codes + cur_ev_codes
    GO_terms = GO_utils.get_go_terms_descendants(biomart_file_path, gene2go_file_path, gene_count_file_path, obo_file_path, ev_codes=ev_codes)
    GO_terms = GO_utils.sort_go_terms(GO_terms)

    term = GO_terms[0]
    utils.predict(term, NUM_FEATURES, rpkm_file_path)
    '''
    for t in GO_terms[0:10]:
        print t.id, ' ', len(t.genes)

    for term in GO_terms:
        utils.predict(term)
        break
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
    '''
    num_folds = 5   # number of folds to use for cross-validation
    loss_function = 'l2'  # Loss function to use. Must be either 'l1' or 'l2'
    logreg_cv_L2 = linear_model.LogisticRegressionCV(cv=num_folds, penalty=loss_function)
    logreg_cv_L2.fit(train.gene_features, train.labels)
    pred_lr_cv_L2 = logreg_cv_L2.predict(test.gene_features)
    utils.print_prediction_results('Cross-Validated Logistic Regression', test.labels, pred_lr_cv_L2,
                              other_info='Norm: ' + loss_function + ', # of Folds: ' + str(num_folds))

    '''

    '''
    svc = SVC(kernel='rbf')
    Cs = np.logspace(-2, 4, 10)
    clf = GridSearchCV(estimator=svc, param_grid=dict(C=Cs),n_jobs=-1)
    clf.fit(train.gene_features, train.labels)
    print 'Best score: ', clf.best_score_
    best_C = clf.best_estimator_.C
    print 'Best C: ', best_C

    clf = SVC(kernel='rbf', C=best_C)
    clf.fit(train.gene_features, train.labels)
    pred_svm = clf.predict(test.gene_features)
    utils.print_prediction_results('SVM', test.labels, pred_svm)
    '''







