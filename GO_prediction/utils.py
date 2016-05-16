
import numpy as np
import random
from math import sqrt, ceil
from sklearn.metrics import roc_auc_score
import GO_utils
from sklearn import linear_model

class GeneData:
    def __init__(self, num_features):
        self.gene_features = np.empty((0, num_features))
        self.labels = []
        self.gene_ids = []

    def append_example(self, example, label, gene_id):
        self.gene_features = np.append(self.gene_features, [example], axis=0)
        self.labels.append(label)
        self.gene_ids.append(gene_id)


def split_data(gene_features, labels, gene_ids_ordered, train_set_size=0.7):
    # train_set_size: Fraction of genes used for training set

    num_examples = len(labels)
    num_features = gene_features.shape[1]
    num_train_examples = int(ceil(train_set_size*num_examples))
    train_indeces = random.sample(range(0, num_examples), num_train_examples)

    train = GeneData(num_features)
    test = GeneData(num_features)

    print 'num examples: ', num_examples
    for idx in range(0, num_examples):
        if idx in train_indeces:
            train.append_example(gene_features[idx], labels[idx], gene_ids_ordered[idx])
        else:
            test.append_example(gene_features[idx], labels[idx], gene_ids_ordered[idx])

    print 'Dimensionality of training set: ', train.gene_features.shape
    print 'Dimensionality of test set: ', test.gene_features.shape
    return train, test


def print_prediction_results(model, fit, labels, predictions, other_info=None):
    print 20*'-'
    print model
    print 20*'-'
    print 'ROC AUC Score: ', roc_auc_score(labels, predictions)
    print fit.coef_
    false_positives = 0
    false_negatives = 0
    for label, pred in zip(labels, predictions):
        if label == 0 and pred == 1:
            false_positives += 1
        elif label == 1 and pred == 0:
            false_negatives += 1

    num_pos_predictions = np.sum(predictions)
    num_neg_predictions = len(predictions) - num_pos_predictions
    print 'False positive rate: ', 1.0 * false_positives / num_pos_predictions
    print 'False negative rate: ', 1.0 * false_negatives / num_neg_predictions

    if other_info:
        print other_info
    print 'done'
    return


def save_prediction_results(fname, GO_id, model, test_data, preds, n_folds=None,
                            tissue_set=None, costs=None, best_cost=None, other_info=None):
    """
    Save results of prediction for a given GO term to a text file.

    :param fname: Name of output file
    :param GO_id: ID of the GO term
    :param model: Type of model used (e.g. "Logistic Regression with L1 Norm")
    :param test_data: GeneData object containing test data
    :param preds: Binary predictions for test data
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
    auc_score = roc_auc_score(test_data.labels, preds)

    print 'auc: ', auc_score

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

    # Writ out the predictions
    out_file.write('# Gene ID\tLabel\tPrediction\n')
    for id,label,pred in zip(test_data.gene_ids, test_data.labels, preds):
        out_file.write(id + '\t' + str(label) + '\t' + str(pred) + '\n')

    out_file.close()



def rand_sample_exclude(li, num_samples, exclude=None):
    """
    Generates a list of unique randomly sampled values from |list|
    Values provided in the list |exclude| will not be selected.

    :param li: The list to sample from
    :param num_samples: Number of samples to take
    :param exclude: List of samples that should not be taken
    :return: The list of |num_samples| samples
    """
    samples = random.sample(li, num_samples)

    if exclude:
        # Check list for any values that should be excluded
        for sample in samples:
            if sample in exclude:
                samples.remove(sample)

        while len(samples) < num_samples:
            sample = random.choice(li)
            if sample not in samples and sample not in exclude:
                samples.append(sample)

    return samples


def logistic_regresssion_L1(term, train, test):
    num_folds = 3   # number of folds to use for cross-validation
    loss_function = 'l1'  # Loss function to use. Must be either 'l1' or 'l2'
    costs = np.logspace(-4, 4, 20)  # 10^(-start) to 10^stop in 10 logarithmic steps
    logreg_cv_L1 = linear_model.LogisticRegressionCV(Cs=costs, cv=num_folds, penalty=loss_function, scoring='roc_auc', solver='liblinear', tol=0.0001)
    logreg_cv_L1.fit(train.gene_features, train.labels)
    best_c = logreg_cv_L1.C_
    pred_lr_cv_L1 = logreg_cv_L1.predict(test.gene_features)
    print_prediction_results('Cross-Validated Logistic Regression', logreg_cv_L1, test.labels, pred_lr_cv_L1,
                                   other_info='Norm: ' + loss_function + ', # of Folds: ' + str(num_folds))
    print 'Best cost (after inverting to obtain true cost): ', 1.0 / best_c

    # Save results
    out_fname = '../data/result_logreg_' + term.id + '.txt'
    save_prediction_results(out_fname, term.id, 'Logistic Regresion with L1 Penalty',
                            test, pred_lr_cv_L1, num_folds,
                            tissue_set=None, costs=costs, best_cost=best_c)


def predict(term, num_features, rpkm_path):
    """
    Predict GO associations for |term|.
    :param term: The GO term
    :return: None
    """

    ensembl_ids = term.genes
    ens_ids_dict = {}
    for id in ensembl_ids:
        ens_ids_dict[id] = True

    print 'Analyzing GO term: ', term.id
    print 'This term has ', len(ensembl_ids), ' genes associated with it.'
    print len(set(ensembl_ids))
    # 1st Pass Through Dataset: Obtain positive training examples
    gene_features, positive_example_rows, gene_ids_ordered, num_transcripts = \
        GO_utils.get_positive_examples(rpkm_path, ens_ids_dict, num_features)

    print 'After pass 1 (inserted positive examples), gene feature matrix has dimension: ', gene_features.shape
    num_positive_examples = len(positive_example_rows)
    num_negative_examples = num_positive_examples
    num_examples = num_positive_examples + num_negative_examples
    print 'num pos: ', num_positive_examples
    print 'num neg: ', num_negative_examples

    # 2nd Pass through dataset: Obtain an equal number of negative training exmaples
    neg_rows = rand_sample_exclude(range(0, num_transcripts), num_negative_examples, exclude=positive_example_rows)

    gene_features_neg, gene_ids_ordered_neg = \
        GO_utils.get_negative_examples(rpkm_path, neg_rows, num_features)
    gene_features = np.append(gene_features, gene_features_neg, axis=0)
    gene_ids_ordered += gene_ids_ordered_neg

    print 'After pass 2 (inserted negative examples), gene feature matrix has dimension: ', gene_features.shape

    # Vector of labels for each example
    labels = num_positive_examples * [1] + num_negative_examples * [0]

    train, test = split_data(gene_features, labels, gene_ids_ordered, train_set_size=0.7)

    # Run Logistic Regression
    logistic_regresssion_L1(term, train, test)





