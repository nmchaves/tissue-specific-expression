"""

    for donor in donors
        for tissue in donor's list:
            write mean exp level in this tissue

"""

import numpy as np

def writeLine(fd, vals, delimiter):
    for val in vals[:-1]:
        matrix_file.write(str(val) + delimiter)
    fd.write(str(vals[-1]) + '\n')

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

    tissues = sorted(getArrayFromFile('../data/tissues.txt'))
    num_tissues = len(tissues)
    tissueToColumn = {}
    for (i, tissue) in enumerate(tissues):
        tissueToColumn[tissue] = i

    # TODO: Also consider sorting matrix by gender or something
    donors = sorted(getArrayFromFile('../data/donors.txt'))

    matrix_file = open('../data/mean_expression_matrix.txt', 'w')
    header = 'Donor_ID\t'
    for tissue in tissues[0:-1]:
        header += tissue + '\t'
    header += tissues[-1] + '\n'
    matrix_file.write(header)

    num_header_lines = 4  # Number of header lines in metadata files
    for donor in donors:
        donor_exp_levels = np.zeros(num_tissues)
        donor_file_path = '../data/Donor_Metadata_Enhanced/donor_meta_' + donor + '.txt'

        with open(donor_file_path) as donor_file:
            for _ in xrange(num_header_lines):
                next(donor_file)
            for line in donor_file:
                vals = line.split('\t')
                tissue = vals[1]
                mean = vals[2]
                col = tissueToColumn[tissue]
                donor_exp_levels[col] = mean

        matrix_file.write(donor + '\t')
        writeLine(matrix_file, donor_exp_levels, '\t')
        donor_file.close()

    matrix_file.close()
