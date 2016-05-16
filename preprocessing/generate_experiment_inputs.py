"""
    This file goes through all GO terms of interest and creates the positive and negative
    examples for running a prediction problem. A positive example is a gene that is known
    to be associated with the GO term, while a negative example is a randomly sampled gene
    from the set of genes that are not known to be associated with the GO term.

"""

from GO_Evidence_Codes import EvidenceCodes
from os import remove

import sys
sys.path.insert(0, '../GO_prediction')
import GO_utils
import utils


def get_all_go_terms():
    ev_codes_obj = EvidenceCodes()
    ev_codes = ev_codes_obj.get_codes(['exp', 'compan', 'auth', 'cur'])
    GO_terms = GO_utils.get_go_terms_descendants(biomart_file_path, gene2go_file_path, gene_count_file_path, obo_file_path, ev_codes=ev_codes)
    GO_terms = GO_utils.sort_go_terms(GO_terms)
    return GO_terms


def get_exp_levels_str(line):

    # Return every character after the 4th tab
    tabs_seen = 0
    for (idx, c) in enumerate(line):
        if c == '\t':
            tabs_seen += 1
            if tabs_seen == 4:
                break

    return line[idx+1:]


def get_ensembl_id(line):
    tab1_index = line.find('\t')
    tab2_index = line.find('\t', tab1_index+1)
    cur_ens_id = line[tab1_index+1:tab2_index]
    # Remove decimal from the ensembl ID
    if '.' in cur_ens_id:
        cur_ens_id = cur_ens_id[0:cur_ens_id.index('.')]

    return cur_ens_id


def make_neg_files(term, num_files, num_examples, num_transcripts, pos_ex_rows):
    """
    Make several different files each containing a randomly sampled set of genes from the set of
    genes in the rpkm file that are not known to be associated with |GOterm|

    :param GOterm: The GO term of interest
    :param num_files: The number of negative files to make.
    :return: None
    """

    # Create Files Containing Negative Examples
    for i in range(0, num_files):

        neg_rows = utils.rand_sample_exclude(range(0, num_transcripts), num_examples, exclude=pos_ex_rows)

        f_name = '../data/experiment_inputs/' + str(term.id) + '_neg_' + str(i) + '.txt'
        neg_file = open(f_name, 'w')

        # Write Headers
        header_1 = '# Negative Examples For GO Term: ' + str(term.id) + '\n'
        header_2 = '# Gene ID' + '\t' + 'Expression Profile\n'
        neg_file.write(header_1 + header_2)

        rpkm_file = open(rpkm_file_path)
        i = 0
        firstLine = True
        for line in rpkm_file:
            if firstLine:
                firstLine = False
                continue

            if i in neg_rows:
                cur_ens_id = get_ensembl_id(line)
                exp_levels_str = get_exp_levels_str(line)
                neg_file.write(cur_ens_id + '\t' + exp_levels_str)

            i += 1

        neg_file.close()
    return


def make_pos_file(term, min_genes, ens_ids_dict):
    """
    Make file containing the genes in the rpkm file that are known to be associated
    with |GOterm|

    :param GOterm: GO_utils.GOterm object containing the term of interest
    :return: True if the file was created, False o/w. In particular, return False
    if there are fewer than |min_genes| in the rpkm file that are associated with
    this GO term.
    """

    f_name = '../data/experiment_inputs/' + str(term.id) + '_pos.txt'
    pos_file = open(f_name, 'w')

    # Write Headers
    header_1 = '# Positive Examples For GO Term: ' + str(term.id) + '\n'
    header_2 = '# Gene ID' + '\t' + 'Expression Profile\n'
    pos_file.write(header_1 + header_2)

    positive_example_rows = []
    genes_added = []

    rpkm_file = open(rpkm_file_path)
    i = 0
    firstLine = True
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            continue

        cur_ens_id = get_ensembl_id(line)

        if cur_ens_id in ens_ids_dict:
            # The IF condition below prevents using the same gene for multiple
            # features. TODO: better method for accounting for multiple transcripts
            # mapping to same gene.
            if cur_ens_id not in genes_added:
                positive_example_rows.append(i)
                genes_added.append(cur_ens_id)
                pos_file.write(cur_ens_id + '\t')
                pos_file.write(get_exp_levels_str(line))
        i += 1

    pos_file.close()
    num_genes_added = len(genes_added)
    print '# of genes actually added from rpkm file: ', num_genes_added
    if num_genes_added >= min_genes:
        return i, positive_example_rows
    else:
        # Delete the file and return None to indicate there weren't enough associations
        remove(f_name)  # TODO: test this
        return None, None


'''

*********************
        Main
*********************
'''
if __name__ == "__main__":

    gene2go_file_path = '../data/gene2go.txt' # If file doesn't exist, then run gene2go = download_ncbi_associations()
    rpkm_file_path = '../../CS341_Data/transcript_rpkm_in_go_nonzero_exp.txt'
    gene_count_file_path = '../data/supp_GO_term_gene_counts.txt'
    biomart_file_path = '../data/biomart_ensembl_to_entrez.txt'
    obo_file_path = '../data/go-basic.obo'

    # TODO: consider doing cases where you use multiple version of positive files instead of
    # just selecting the first gene each time.
    num_negative_files = 1  # TODO: use more negative files

    # A given GO term must be associated with at least |min_gene_associations| number of
    # genes that are present in the rpkm file for it to be used in prediction.
    min_gene_associations = 10
    num_usable_terms = 0

    GO_terms = get_all_go_terms()
    for (t, term) in enumerate(GO_terms):

        print t, 'th term. ID: ', term.id, ' # of genes (before going through rpkm): ', len(term.genes)

        if len(term.genes) < min_gene_associations:
            break

        ensembl_ids = term.genes
        ensembl_ids_dict = {}
        for id in ensembl_ids:
            ensembl_ids_dict[id] = True

        num_transcripts, pos_rows = make_pos_file(term, min_gene_associations, ensembl_ids_dict)
        if pos_rows is None:
            # make_pos_file returns None when there were fewer than |min_gene_associations|
            # so this GO term will not be used
            continue

        make_neg_files(term, num_negative_files, len(pos_rows), num_transcripts, pos_rows)

        num_usable_terms += 1
        print '# of usable GO terms found so far: ', num_usable_terms
