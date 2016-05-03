'''
    Consider biological process X
    Training example: gene
    Label (binary): 1 if gene is associated with X in GO,
    0 otherwise. Can obtain negative examples by randomly sampling
    the genes that are not known to be associated with X

    # TODO: cross-validation, AUC score, clean up pipeline, add in tissue
    selection to the pipeline

    # TODO: use only 53 tissues
'''

import numpy as np
import random
from sklearn import linear_model
from goatools.associations import read_ncbi_gene2go
from goatools.base import download_ncbi_associations
import math
from sklearn.metrics import roc_auc_score
from sklearn.metrics import mean_squared_error
from sklearn.svm import SVC

#from goatools.base import download
#obo_fname = download_go_basic_obo()

def rand_sample_exclude(list, num_samples, exclude=None):
    '''
    Generates a list of unique randomly sampled values from |list|

    :param list: The list to sample from
    :param num_samples: Number of samples to take
    :param exclude: List of samples that should not be taken
    :return: The list of |num_samples| samples
    '''
    samples = random.sample(list, num_samples)

    if exclude:
        # Check list for any values that should be excluded
        for sample in samples:
            if sample in exclude:
                samples.remove(sample)

        while len(samples) < num_samples:
            sample = random.choice(list)
            if sample not in samples and sample not in exclude:
                samples.append(sample)

    return samples

def map_entrez_to_ensembl(path):
    dict = {}
    file = open(path)
    for line in file:
        vals = line.split('\t')
        ens_gene_id = vals[0]
        entrez_id = vals[2]
        dict[entrez_id] = ens_gene_id

    file.close()
    return dict


def get_ensembl_ids(go_process_id, biomart_fpath, ev_codes=None):

    entrez_to_ensembl = map_entrez_to_ensembl(biomart_fpath)

    gene2go = 'data/gene2go.txt' # If file doesn't exist, then replace this line with gene2go = download_ncbi_associations()

    # taxids=[9606] means select only human.
    go_to_entrez_ids_human = read_ncbi_gene2go(gene2go, taxids=[9606], go2geneids=True)
    print("{N} GO terms associated with human NCBI Entrez GeneIDs".format(N=len(go_to_entrez_ids_human)))

    entrez_ids = go_to_entrez_ids_human[GO_PROCESS_ID]
    print '# of Entrez IDs associated with ', GO_PROCESS_ID, ' = ', len(entrez_ids)
    ensembl_ids = []
    for ent_id in entrez_ids:
        if str(ent_id) in entrez_to_ensembl:
            ensembl_ids.append(entrez_to_ensembl[str(ent_id)])

    print '# of Ensembl IDs associated with ', GO_PROCESS_ID, ' = ', len(ensembl_ids)
    return ensembl_ids


def get_positive_examples(rpkm_path, ens_ids):

    gene_features = np.empty((0, NUM_SAMPLES))
    positive_example_rows = []
    gene_ids_ordered = []
    i = 0
    rpkm_file = open(rpkm_file_path)
    firstLine = True
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            continue

        tab1_index = line.find('\t')
        tab2_index = line.find('\t', tab1_index+1)
        cur_ens_id = line[tab1_index+1:tab2_index]
        # Remove decimal from the ensembl ID
        if '.' in cur_ens_id:
            cur_ens_id = cur_ens_id[0:cur_ens_id.index('.')]

        if cur_ens_id in ens_ids:
            # This is IF condition prevents using the same gene for multiple
            # features. TODO: better method for accounting for multiple transcripts
            # mapping to same gene.
            if cur_ens_id not in gene_ids_ordered:
                positive_example_rows.append(i)
                gene_ids_ordered.append(cur_ens_id)
                exp_levels_str = line.rstrip().split('\t')[4:]
                exp_levels = [float(exp_level) for exp_level in exp_levels_str]
                gene_features = np.append(gene_features, [exp_levels], axis=0)
        i += 1
    rpkm_file.close()

    return gene_features, positive_example_rows, gene_ids_ordered, i


def get_negative_examples(rpkm_path, neg_ex_rows):

    gene_features_neg = np.empty((0, NUM_SAMPLES))
    gene_ids_ordered_neg = []
    rpkm_file = open(rpkm_file_path)
    i = 0
    firstLine = True
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            continue

        if i in neg_ex_rows:
            vals = line.rstrip().split('\t')
            cur_ens_id = vals[1]

            # Remove decimal from the ensembl ID
            if '.' in cur_ens_id:
                cur_ens_id = cur_ens_id[0:cur_ens_id.index('.')]
            gene_ids_ordered_neg.append(cur_ens_id)
            exp_levels_str = vals[4:]
            exp_levels = [float(exp_level) for exp_level in exp_levels_str]
            gene_features_neg = np.append(gene_features_neg, [exp_levels], axis=0)
        i += 1

    rpkm_file.close()
    return gene_features_neg, gene_ids_ordered_neg

def print_prediction_results(model, labels, predictions, other_info=None):
    print 20*'-'
    print model
    print 20*'-'
    print 'Root Mean Square Error: ', math.sqrt(mean_squared_error(labels, predictions))
    print 'ROC AUC Score: ', roc_auc_score(labels, predictions)
    false_positives = 0
    false_negatives = 0
    for label, pred in zip(labels, predictions):
        if label == 0 and pred == 1:
            false_positives += 1
        elif label == 1 and pred == 0:
            false_negatives += 1

    print 'False positive rate: ', 1.0 * false_positives / len(labels)
    print 'False negative rate: ', 1.0 * false_negatives / len(labels)

    if other_info:
        print other_info


"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    GO_PROCESS_ID = 'GO:0001889'  # Biological Process ID in Gene Ontology
    biomart_file_path = 'data/biomart_ensembl_to_entrez.txt'
    #rpkm_file_path = '../../../Documents/Stanford/CS341_Data/transcript_rpkm_top_10000_var.txt'
    rpkm_file_path = '../../../Downloads/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt'
    go_evidence_codes = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']
    ensembl_ids = get_ensembl_ids(GO_PROCESS_ID, biomart_file_path, ev_codes=go_evidence_codes)

    NUM_SAMPLES = 8555

    # 1st Pass Through Dataset: Obtain positive training examples
    gene_features, positive_example_rows, gene_ids_ordered, num_transcripts = get_positive_examples(rpkm_file_path, ensembl_ids)

    print 'After pass 1 (inserting positive examples), gene feature matrix has dimension: ', gene_features.shape
    num_positive_examples = len(positive_example_rows)
    num_negative_examples = num_positive_examples
    num_examples = num_positive_examples + num_negative_examples

    # 2nd Pass through dataset: Obtain an equal number of negative training exmaples
    max_row_id = num_transcripts-1
    negative_example_rows = rand_sample_exclude(range(0, max_row_id+1), num_positive_examples, exclude=positive_example_rows)

    gene_features_neg, gene_ids_ordered_neg = get_negative_examples(rpkm_file_path, negative_example_rows)
    gene_features = np.append(gene_features, gene_features_neg, axis=0)
    gene_ids_ordered += gene_ids_ordered_neg

    print 'After pass 2 (inserting negative examples), gene feature matrix has dimension: ', gene_features.shape

    # Vector of labels for each example
    labels = num_positive_examples * [1] + num_negative_examples * [0]

    # Split into training and test sets
    # gene_features_train, labels_train, gene_features_test, labels_test =
    # split_to_train_and_test(gene_features, TRAIN_SET_SIZE)
    TRAIN_SET_SIZE = 0.7  # Fraction of genes used for training set
    num_train_examples = int(math.ceil(TRAIN_SET_SIZE*num_examples))
    train_indeces = random.sample(range(0, num_examples), num_train_examples)

    gene_features_train = np.empty((0, NUM_SAMPLES))
    gene_features_test = np.empty((0, NUM_SAMPLES))
    labels_train = []
    labels_test = []
    gene_ids_ordered_train = []
    gene_ids_ordered_test = []
    print 'num ex: ', num_examples
    for idx in range(0, num_examples):
        if idx in train_indeces:
            gene_features_train = np.append(gene_features_train, [gene_features[idx]], axis=0)
            labels_train.append(labels[idx])
            gene_ids_ordered_train.append(gene_ids_ordered[idx])
        else:
            gene_features_test = np.append(gene_features_test, [gene_features[idx]], axis=0)
            labels_test.append(labels[idx])
            gene_ids_ordered_test.append(gene_ids_ordered[idx])

    print 'Dimensionality of training data: ',gene_features_train.shape

    logreg = linear_model.LogisticRegression(C=1e5)
    logreg.fit(gene_features_train, labels_train)

    pred = logreg.predict(gene_features_test)
    print_prediction_results('Logistic Regression', labels_test, pred)

    # Logistic Regression with Cross-Validation
    #costs = []
    num_folds = 10   # number of folds to use for cross-validation
    loss_function = 'l1'  # Loss function to use. Must be either 'l1' or 'l2'
    logreg_cv = linear_model.LogisticRegressionCV(cv=num_folds, penalty=loss_function)
    logreg.fit(gene_features_train, labels_train)

    pred = logreg.predict(gene_features_test)
    print_prediction_results('Cross-Validated Logistic Regression', labels_test, pred,
                              other_info='Norm: ' + loss_function + ', # of Folds: ' + str(num_folds))

    num_folds = 10   # number of folds to use for cross-validation
    loss_function = 'l2'  # Loss function to use. Must be either 'l1' or 'l2'
    logreg_cv = linear_model.LogisticRegressionCV(cv=num_folds, penalty=loss_function)
    logreg.fit(gene_features_train, labels_train)

    pred = logreg.predict(gene_features_test)
    print_prediction_results('Cross-Validated Logistic Regression', labels_test, pred
                             , other_info='Norm: ' + loss_function + ', # of Folds: ' + str(num_folds))

    # SVM
    clf = SVC()
    clf.fit(gene_features_train, labels_train)

    pred_svm = clf.predict(gene_features_test)
    print_prediction_results('SVM', labels_test, pred)