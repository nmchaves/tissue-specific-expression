import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile, join
from sklearn.metrics import roc_auc_score, roc_curve, auc
import pandas as pd
from collections import defaultdict
from scipy.stats import cumfreq


def get_GO_gene_counts(input_file):
    counts_dict = {}
    counts_file = open(input_file)
    for (i, line) in enumerate(counts_file):
        if i < 2:
            continue
        data = line.rstrip().split('\t')
        counts_dict[data[0]] = int(data[1])
    counts_file.close()
    return counts_dict


def get_prediction_results(results_dir):
    results_files = [f for f in listdir(results_dir) if isfile(join(results_dir, f))]

    GO_terms = {}
    roc_auc_scores = []
    roc_auc_score_line = 2
    for rf in results_files:
        f = open(results_dir + '/' + rf)
        GO_term = None
        labels = []
        preds = []
        dec_func_scores = []
        probs = []
        for (i, line) in enumerate(f):
            if i == 0:
                GO_term = line.rstrip().split(' ')[-1]
            elif i == roc_auc_score_line:
                vals = line.rstrip().split(' ')
                roc_auc_scores.append(float(vals[-1]))
            elif line[0:4] == 'ENSG':
                vals = line.rstrip().split('\t')
                labels.append(int(vals[1]))
                preds.append(int(vals[2]))
                dec_func_scores.append(float(vals[3]))
                probs.append(float(vals[4]))

        GO_terms[GO_term] = (labels, preds, dec_func_scores, probs)
    return GO_terms, roc_auc_scores

def get_GO_gene_cnts():
    cnts_file = open('../data/GO_terms_final_gene_counts.txt')
    GO_to_num_genes = {}
    for (i, line) in enumerate(cnts_file):
        if i < 2:
            continue
        vals = line.rstrip().split('\t')
        GO_id = vals[0]
        num_genes = int(vals[1])
        GO_to_num_genes[GO_id] = num_genes
    return GO_to_num_genes


def plot_roc(fprs, tprs, title):
    plt.figure()
    for (fpr, tpr) in zip(fprs, tprs):
        plt.plot(fpr, tpr, 'gray')
    plt.plot([0, 1], [0, 1], 'k--')  # Plot the 50% line
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc="lower right")
    plt.show()
    return


def plot_roc_heat(fprs, tprs, gene_cnts, title):

    NCURVES = len(gene_cnts)
    xs = fprs

    fig = plt.figure()
    ax = fig.add_subplot(111)

    cmap = plt.cm.YlOrRd
    min_count = int(min(gene_cnts))
    max_count = int(max(gene_cnts))
    print min_count, max_count
    scalarMap = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min_count, vmax=max_count))

    for idx in range(NCURVES):
        colorVal = scalarMap.to_rgba(gene_cnts[idx])
        ax.plot(fprs[idx], tprs[idx], color=colorVal)

    ax.plot([0, 1], [0, 1], 'black')  # Plot the 50% line

    # Generate colorbar
    scalarMap.set_array([])  # You have to set a dummy-array for this to work...
    cbar = plt.colorbar(scalarMap)
    cbar.set_label('# of Genes')
    #cbar.set_ticks(gene_cnts)
    #cbar.set_ticklabels(['{:4.1f}'.format(yi) for yi in y]) # Make 'em nicer-looking

    ax.grid()
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.show()
    return


def make_roc_curves(GO_terms_map, GO_cnts):
    false_pos_rates = []
    true_pos_rates = []
    GO_terms_list = GO_terms_map.keys()
    gene_counts = [GO_cnts[term] for term in GO_terms_list]
    sorted_tuples = sorted(zip(GO_terms_list, gene_counts), key=lambda tup: tup[1], reverse=False)
    GO_terms_list = [tup[0] for tup in sorted_tuples]
    gene_counts = [tup[1] for tup in sorted_tuples]
    auc_scores = []
    for term in GO_terms_list:
        (cur_labels, cur_preds, cur_decs, cur_probs) = GO_terms_map[term]
        cur_fpr, cur_tpr, _ = roc_curve(cur_labels, cur_probs)

        #cur_fpr, cur_tpr, _ = roc_curve(cur_labels, cur_decs)
        false_pos_rates.append(cur_fpr)
        true_pos_rates.append(cur_tpr)
        auc_scores.append(roc_auc_score(cur_labels, cur_probs))

    plot_roc(false_pos_rates, true_pos_rates, 'ROC for Log Transformed Expressions')
    plot_roc_heat(false_pos_rates, true_pos_rates, gene_counts, 'ROC for Log Transformed Expressions')

    plt.plot(gene_counts,auc_scores, "o")
    plt.xlabel('Gene Counts')
    plt.ylabel('AUC Score')
    ax = plt.gca()
    ax.grid()
    plt.show()

    return auc_scores


def get_go_terms():
    f_name = '../data/GO_terms_final_gene_counts.txt'
    GO_counts_file = open(f_name)

    terms = []
    for (i, line) in enumerate(GO_counts_file):
        if i < 2:
            continue
        term = line.split('\t')[0]
        terms.append(term)
    return terms

def get_tissue_list(tissue_fpath):
    tissue_file = open(tissue_fpath)
    for line in tissue_file:
        tissues = line.rstrip().split('\t')
        break
    return tissues

def get_1_tissue_aucs(GO_term, tissue_list, loss=None):
    """
    This function gets the AUC scores of predicting the gene associations of
    |GO_term| where each prediction task only used features from an individual
    tissue.

    :param GO_term: The GO term
    :param tissue_list: List of tissues that were used for separate prediction
    tasks. If len(tissue_list)=53, then we performed 53 separate prediciton
    problems, where each prediction problem used a different 1 of the 53 tissues.
    :return: A list of AUC scores in the same order as |tissue_list|. The ith element in this
    list is the AUC score for predicting this GO term using only tissue i.
    """

    aucs_1_tissue = []
    if loss:
        results_dir = 'pca_results_1_tissue_' + loss + '/' + GO_term + '/'
    else:
        results_dir = 'results_1_tissue/' + GO_term + '/'
    for tissue in tissue_list:
        # Get the AUC score when using features from only this tissue
        f_name = results_dir + 'logreg_' + tissue + '.txt'
        rf = open(f_name)
        for (i, line) in enumerate(rf):
            if i == 2:
                auc_score = float(line.split(' ')[-1])
                aucs_1_tissue.append(auc_score)
            elif i > 2:
                break
    return aucs_1_tissue

def get_all_1_tissue_aucs(GO_terms, tissue_list, loss=None):
    aucs = np.zeros(shape=(len(GO_terms),len(tissue_list)))  # aucs[i][j] is auc score for using jth tissue features to predict ith GO term
    for (i, term) in enumerate(GO_terms):
        aucs[i, :] = get_1_tissue_aucs(term, tissue_list, loss)
    return aucs


def map_GO_to_GTEX():
    inputFilename = '../data/GO_terms_final_gene_counts.txt'
    GO_list_file = open(inputFilename)
    GO_list = np.loadtxt(GO_list_file,skiprows=2,usecols=[0],dtype='S10',delimiter='\t')

    inputFilename = '../data/Tissue_Name_Mappings.csv'
    tissue_data = pd.read_csv(inputFilename,header=None)
    map_BTO_to_GTEX = defaultdict(list)

    for index,row in tissue_data.iterrows():
        GTEX_tissue = row[0]
        BTO_tissues = row[1:]
        for tissue in BTO_tissues.dropna():
            map_BTO_to_GTEX[tissue].append(GTEX_tissue)

    inputFilename = '../data/BTO_GO.csv'
    BTO_data = pd.read_csv(inputFilename,skiprows=[0])
    map_GO_to_GTEX = defaultdict(list)

    for index,row in BTO_data.iterrows():
        tissue = row[1]
        if tissue in map_BTO_to_GTEX:
            GO_IDs = row[2:]
            for GO_ID in GO_IDs.dropna():
                if GO_ID in GO_list:
                    map_GO_to_GTEX[GO_ID] = list(set(map_GO_to_GTEX[GO_ID] + map_BTO_to_GTEX[tissue]))

    #inputFile.close()
    return map_GO_to_GTEX

def map_GTEX_to_GO(map_GO_to_GTEX):
    GTEX_to_GO = defaultdict(list)
    for GO_ID,tissues in map_GO_to_GTEX.items():
        for tissue in tissues:
            if GO_ID not in GTEX_to_GO[tissue]:
                GTEX_to_GO[tissue].append(GO_ID)
    return GTEX_to_GO


def map_GTEX_to_cols(dir_path, tissues):
    GTEX_to_samples = {}
    for tissue in tissues:
        cols = []
        fpath = dir_path + 'tissue_meta_' + tissue + '.txt'
        meta_file = open(fpath)
        for (i, line) in enumerate(meta_file):
            if i < 1:
                continue
            else:
                cols.append(int(line.split('\t')[0]))
        GTEX_to_samples[tissue] = cols
        meta_file.close()
    return GTEX_to_samples


def make_boxplot(vals, labels):
    plt.figure(figsize=(18, 6))
    plt.margins(0.01)
    plt.ylim([0, 1])
    ax = plt.gca()
    ax.xaxis.grid(which='both')
    plt.xticks(range(len(labels)), labels, rotation='vertical')
    ax.boxplot(vals, labels=labels)
    plt.show()
    return



