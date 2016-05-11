
import numpy as np
import random
from math import sqrt, ceil
from sklearn.metrics import roc_auc_score


class GeneData:
    def __init__(self, num_features):
        self.gene_features = np.empty((0, num_features))
        self.labels = []
        self.gene_ids_ordered = []

    def append_example(self, example, label, gene_id):
        self.gene_features = np.append(self.gene_features, [example], axis=0)
        self.labels.append(label)
        self.gene_ids_ordered.append(gene_id)


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


def print_prediction_results(model, labels, predictions, other_info=None):
    print 20*'-'
    print model
    print 20*'-'
    print 'ROC AUC Score: ', roc_auc_score(labels, predictions)
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
                            tissue_set=None, best_cost=None, other_info=None):
    """
    Save results of prediction for a given GO term to a text file.

    :param fname: Name of output file
    :param GO_id: ID of the GO term
    :param model: Type of model used (e.g. "Logistic Regression with L1 Norm")
    :param test_data: GeneData object containing test data
    :param preds: Binary predictions for test data
    :param n_folds: # of folds used for cross-validation
    :param tissue_set: Set of tissues used. If no option is specified, it's assumed that all tissues were used
    :param best_cost: If applicable, the best cost parameter (as determined by cross-validation)
    :param other_info: Other info to store in header of output file
    :return: None
    """
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
    if best_cost:
        out_file.write('# Best cost parameter (determined by CV): ' + str(best_cost) + '\n')
    if other_info:
        out_file.write('# ' + other_info + '\n')

    # Writ out the predictions
    out_file.write('# Gene ID\tLabel\tPrediction\n')
    for id,label,pred in zip(test_data.gene_ids_ordered, test_data.labels, preds):
        out_file.write(id + '\t' + str(label) + '\t' + str(pred) + '\n')

    out_file.close()
    """
    print 'hi'


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
