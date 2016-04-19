"""
    This file checks that after we generated donor files from
    "GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt", we
    didn't eliminate any donors from the dataset.

    for file in donor_file:
        calc_mean

    sort means

"""


def getArrayFromFile(path):
    """
    This function takes a path to a file containing a tab-delimited
    text file and converts that file into an array.

    :param path: Path to the file
    :return: The text file
    """
    f = open(path)
    for line in f:
        return line.split('\t')

"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    # If you run this yourself, you should set these paths properly!!
    donor_file_path = 'data/donors.txt'
    donor_data_path = '~/Documents/Stanford/CS341_Data/donor_matrices_fixed'

    donors = getArrayFromFile(donor_file_path)
    donor_file_prefix = 'donor_'

    mean_expression_levels = []  # Array containing mean expression level for each user
    i = 0
    for donor in donors:
        file_name = donor_file_prefix + donor + '.txt'
        donor_file = open(file_name)
        mean_expression = i
        i += 1
        mean_expression_levels.append((donor, mean_expression))
        donor_file.close()

    # Sort by expression level
