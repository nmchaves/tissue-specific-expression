"""
    This file performs step 2/3 of preprocessing of the GTEx RPKM file:
    "GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt"

    This should be ran AFTER filter_by_go.py

    This script takes the top |NUM_TRANSCRIPTS_TO_RETAIN| rows by variance and
    saves them into a text file. In other words, this file removes the low variance
    rows from the rpkm file's matrix.

"""

import numpy as np

"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    NUM_TRANSCRIPTS_TO_RETAIN = 10000

    path_to_rpkm_file = '../../../Documents/Stanford/CS341_Data/transcript_rpkm_in_go.txt'

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

    file_name = '../../../Documents/Stanford/CS341_Data/transcript_rpkm_top_' + str(NUM_TRANSCRIPTS_TO_RETAIN) + '_var.txt'
    rpkm_file_top_var = open(file_name, 'w')
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





