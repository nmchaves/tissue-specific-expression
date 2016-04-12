"""
    This file performs all filtering steps of preprocessing of the GTEx RPKM file:
    "GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt".

    It goes through each line of the above rpkm file. If the current transcript
    corresponds to a gene in Gene Ontology, then the current line is written out
    to the file "transcript_rpkm_in_go.txt". Otherwise, the current line is written
    out to the file "transcript_rpkm_not_in_go.txt".

"""

import numpy as np

def filter_by_list(path_to_rpkm_file, path_to_output_file, path_to_list):
    """
    It goes through each line of the above rpkm file. If the current transcript
    corresponds to a gene in Gene Ontology, then the current line is written out
    to the file "transcript_rpkm_in_go.txt". Otherwise, the current line is written
    out to the file "transcript_rpkm_not_in_go.txt".

    :param path_to_rpkm_file: input file name
    :param path_to_output_file: ouput file name
    :param path_to_list: gene conversion file name
    :return: No return value. 
    """
    # path_to_rpkm_file = '../../../Downloads/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt'
    # path_to_go_file   = '../../../Downloads/trans_gene_name_filtered_by_GO.txt'
    # path_to_in_go_file = 'transcript_rpkm_in_go.txt'
    # path_to_not_in_go_file ='transcript_rpkm_not_in_go.txt' 

    # Read in GO file and generate dictionary containing transcripts
    go_file = open(path_to_list)
    transcriptIdsInGO = {}
    firstLine = True
    for line in go_file:
        if firstLine:
            firstLine = False
            continue
        transcriptId = line[0:line.index('\t')]
        transcriptIdsInGO[transcriptId] = True
    go_file.close()

    rpkm_file = open(path_to_rpkm_file)
    rpkm_file_retained = open(path_to_output_file, 'w')
    # rpkm_file_removed = open(path_to_not_in_go_file, 'w')
    firstLine = True
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            # Write the header row into both files
            rpkm_file_retained.write(line)
            # rpkm_file_removed.write(line)
            continue
        else:
            transcriptId = line[0:line.index('\t')]
            if '.' in transcriptId:
                transcriptId = transcriptId[0:transcriptId.index('.')]

            if transcriptId in transcriptIdsInGO:
                rpkm_file_retained.write(line)
            # else:
            #     rpkm_file_removed.write(line)

    rpkm_file_retained.close()
    # rpkm_file_removed.close()


def filter_by_var(path_to_rpkm_file, path_to_output_file, NUM_TRANSCRIPTS_TO_RETAIN):
    """
        This file performs the 2nd step of preprocessing the GTEx RPKM file:
        "GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt"

        This should be ran AFTER filter_by_go.py

        This script takes the top |NUM_TRANSCRIPTS_TO_RETAIN| rows by variance and
        saves them into a text file. In other words, this file removes the low variance
        rows from the rpkm file's matrix.

        :param path_to_rpkm_file: input file name
        :param path_to_output_file: ouput file name
        :param NUM_TRANSCRIPTS_TO_RETAIN: number of transcripts to retain 
        :return: No return value. 

    """

    # NUM_TRANSCRIPTS_TO_RETAIN = 10000
    # Path to the file generated using filter_by_go.py
    # path_to_rpkm_file = '../../../Documents/Stanford/CS341_Data/transcript_rpkm_in_go.txt'
    # path_to_output_file= '../../../Documents/Stanford/CS341_Data/transcript_rpkm_top_' + str(NUM_TRANSCRIPTS_TO_RETAIN) + '_var.txt'

    # First pass: compute variances and record indices sorted by variance
    rpkm_file = open(path_to_rpkm_file)

    variances = []
    i = 0
    firstLine = True
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            continue
        expressionLevels = np.array(line.split('\t')[4:]).astype(np.float)
        variances.append((i, np.var(expressionLevels)))
        i += 1

    # Sort the list by variance (the 2nd element in the tuple)
    variances = sorted(variances, key=lambda tup: tup[1], reverse=True)
    topVarianceIndices = [index for (index, var) in variances[0:NUM_TRANSCRIPTS_TO_RETAIN]]
    topVarianceIndices = sorted(topVarianceIndices)  # Sort by index

    rpkm_file.close()

    # Second pass: write the processed data out to new .rpkm file
    rpkm_file = open(path_to_rpkm_file)

    rpkm_file_top_var = open(path_to_output_file, 'w')
    firstLine = True
    nextIndex = topVarianceIndices.pop(0)
    i = 0
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            rpkm_file_top_var.write(line)
            continue
        else:
            if i == nextIndex:
                rpkm_file_top_var.write(line)
                if len(topVarianceIndices) > 0:
                    nextIndex = topVarianceIndices.pop(0)
                else:
                    break
            i += 1

    rpkm_file.close()
    rpkm_file_top_var.close()


"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    data_dir = '../data/'
    NUM_TRANSCRIPTS_TO_RETAIN = 10000;

    raw_input  = data_dir + 'GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt'
    geneListGO = data_dir + 'trans_gene_name_filtered_by_GO.txt'
    output_go  = data_dir + 'transcript_rpkm_in_go.txt'
    output_var = data_dir + 'transcript_rpkm_top_' + str(NUM_TRANSCRIPTS_TO_RETAIN) + '_var.txt' 

    print 'Filtering by GO gene list'
    filter_by_list(raw_input, output_go, geneListGO)
    print 'Filtering by variance'
    filter_by_var(output_go, output_var, NUM_TRANSCRIPTS_TO_RETAIN)
