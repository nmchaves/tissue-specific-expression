"""
    This file goes through all GO terms of interest and creates the positive and negative
    examples for running a prediction task.

    A positive example is a gene that is known to be associated with the GO term, while
    a negative example is a randomly sampled gene from the set of genes that are not known
    to be associated with the GO term.

"""

from os import remove

from GO_Evidence_Codes import EvidenceCodes
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
        header_2 = '# Row_Idx\tGene_ID\tExpression_Profile\n'
        neg_file.write(header_1 + header_2)

        rpkm_file = open(rpkm_file_path)
        for (i, line) in enumerate(rpkm_file):
            if i == 0:
                continue

            if i-1 in neg_rows:
                cur_ens_id = get_ensembl_id(line)
                exp_levels_str = get_exp_levels_str(line)
                neg_file.write(str(i-1) + '\t' + cur_ens_id + '\t' + exp_levels_str)

        neg_file.close()
    return

def make_neg_file_from(neg_fname, out_fname, term, rpkm):

    neg_rows = []
    neg_file = open(neg_fname)
    for (i, line) in enumerate(neg_file):
        if i < 2:
            continue
        row_idx = int(line.split('\t')[0])
        neg_rows.append(row_idx)
    neg_file.close()

    neg_rows = sorted(neg_rows)

    neg_out = open(out_fname, 'w')
    header_1 = '# Positive Examples For GO Term: ' + term + '\n'
    header_2 = '# Row_Idx\tGene_ID\tExpression_Profile\n'
    neg_out.write(header_1 + header_2)

    rpkm_file = open(rpkm)
    for (i, line) in enumerate(rpkm_file):
        if i == 0:
            continue
        elif i-1 == neg_rows[0]:
            neg_rows.pop(0)
            neg_out.write(str(i-1) + '\t')

            cur_ens_id = get_ensembl_id(line)
            neg_out.write(cur_ens_id + '\t')
            neg_out.write(get_exp_levels_str(line))
            if len(neg_rows) == 0:
                break
    neg_file.close()


def make_pos_file_from(pos_fname, out_fname, term, rpkm):

    pos_rows = []
    pos_file = open(pos_fname)
    for (i, line) in enumerate(pos_file):
        if i < 2:
            continue
        row_idx = int(line.split('\t')[0])
        pos_rows.append(row_idx)
    pos_file.close()

    pos_rows = sorted(pos_rows)

    pos_out = open(out_fname, 'w')
    header_1 = '# Positive Examples For GO Term: ' + term + '\n'
    header_2 = '# Row_Idx\tGene_ID\tExpression_Profile\n'
    pos_out.write(header_1 + header_2)

    rpkm_file = open(rpkm)
    for (i, line) in enumerate(rpkm_file):
        if i == 0:
            continue
        elif i-1 == pos_rows[0]:
            pos_rows.pop(0)
            pos_out.write(str(i-1) + '\t')
            cur_ens_id = get_ensembl_id(line)
            pos_out.write(cur_ens_id + '\t')
            pos_out.write(get_exp_levels_str(line))
            if len(pos_rows) == 0:
                break
    pos_file.close()


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
    header_2 = '# Row_Idx\tGene_ID\tExpression_Profile\n'
    pos_file.write(header_1 + header_2)

    positive_example_rows = []
    genes_added = []

    rpkm_file = open(rpkm_file_path)
    for (i, line) in enumerate(rpkm_file):
        if i == 0:
            continue

        cur_ens_id = get_ensembl_id(line)

        if cur_ens_id in ens_ids_dict:
            # The IF condition below prevents using the same gene for multiple
            # features. TODO: better method for accounting for multiple transcripts
            # mapping to same gene.
            if cur_ens_id not in genes_added:
                positive_example_rows.append(i-1)
                genes_added.append(cur_ens_id)
                pos_file.write(str(i-1) + '\t')
                pos_file.write(cur_ens_id + '\t')
                pos_file.write(get_exp_levels_str(line))

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
    gene_count_file_path = '../data/GO_terms_final_gene_counts.txt'
    biomart_file_path = '../data/biomart_ensembl_to_entrez.txt'
    obo_file_path = '../../CS341_Data/go-basic.obo.txt'

    # TODO: consider doing cases where you use multiple version of positive files instead of
    # just selecting the first gene each time.
    num_negative_files = 1  # TODO: use more negative files

    # A given GO term must be associated with at least |min_gene_associations| number of
    # genes that are present in the rpkm file for it to be used in prediction.
    min_gene_associations = 10
    num_usable_terms = 0

    existing_dir = '../data/experiment_inputs/'

    new_dir = '../data/pca_experiment_inputs/'
    new_rpkm = '../../CS341_Data/log_norm_pca_transcript_rpkm_in_go_nonzero_exp.txt'

    new_dir2 = '../data/median_experiment_inputs/'
    new_rpkm2 = '../../CS341_Data/log_norm_median_transcript_rpkm_in_go_nonzero_exp.txt'

    GO_terms = get_all_go_terms()
    for (t, term) in enumerate(GO_terms):
        print t, 'th term. ID: ', term.id, ' # of genes (before going through rpkm): ', len(term.genes)
        term_id = term.id

        if existing_dir:
            # Generate matching positive/negative sets from the existing ones
            pos_fname = existing_dir + term_id + '_pos.txt'
            new_pos_fname = new_dir + term_id + '_pos.txt'
            make_pos_file_from(pos_fname, new_pos_fname, term_id, new_rpkm)
            neg_fname = existing_dir + term_id + '_neg_0.txt'
            new_neg_fname = new_dir + term_id + '_neg_0.txt'
            make_neg_file_from(neg_fname, new_neg_fname, term_id, new_rpkm)

            new_pos_fname2 = new_dir2 + term_id + '_pos.txt'
            make_pos_file_from(pos_fname, new_pos_fname2, term_id, new_rpkm2)
            new_neg_fname2 = new_dir2 + term_id + '_neg_0.txt'
            make_neg_file_from(neg_fname, new_neg_fname2, term_id, new_rpkm2)


            #out_fname = '../data/pca_experiment_inputs/' + str(term.id) + '_pos.txt'
        else:
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
