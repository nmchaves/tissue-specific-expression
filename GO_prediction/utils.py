
import numpy as np
import random
from math import sqrt, ceil
from sklearn.metrics import roc_auc_score
from sklearn.metrics import mean_squared_error


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
    # Fraction of genes used for training set

    num_examples = len(labels)
    num_features = gene_features.shape[1]
    num_train_examples = int(ceil(train_set_size*num_examples))
    train_indeces = random.sample(range(0, num_examples), num_train_examples)

    train = GeneData(num_features)
    test = GeneData(num_features)

    num_examples = 315

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
    print 'Root Mean Square Error: ', sqrt(mean_squared_error(labels, predictions))
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
