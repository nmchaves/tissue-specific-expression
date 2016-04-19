"""
    This script generates a file containing a tab-delimited matrix
    for each donor.

    The matrix files are generated incrementally. If this program is interrupted at some
    point, then we'll have saved the files for each donor that has been processed so far.

"""

def buildDonorsToColumnsDict(path_to_rpkm_file):
    """
    This function creates a dictionary that maps a donor to the list
    of columns that correspond to this donor in the RPKM file.

    :param path_to_rpkm_file: Path to the RPKM txt file
    :return: The dictionary
    """

    rpkm_file = open(path_to_rpkm_file)

    donorsDict = {}
    for line in rpkm_file:
        lineAsArr = line.split('\t')
        for (col, sampleId) in enumerate(lineAsArr):
            if col >= 4:
                donorId = sampleId.split('-')[1]
                if donorId not in donorsDict:
                    donorsDict[donorId] = [col]
                else:
                    donorsDict[donorId].append(col)
        break  # Only examine the 1st line of rpkm file

    rpkm_file.close()
    return donorsDict


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

    donors = getArrayFromFile('../data/donors.txt')
    rpkm_file_path = '../../../../Documents/Stanford/CS341_Data/transcript_rpkm_top_10000_var.txt'
    donorsToColumns = buildDonorsToColumnsDict(rpkm_file_path)

    def generate_donor_matrix(donorId):
        """
        This function generates the donor matrix for |donorId|.
        After this function completes there will be a file containing
        a tab-delimited matrix. The matrix columns correspond to samples and the
        matrix rows correspond to target IDs. The cell values are expression
        levels.

        :param donorId: The donor's ID
        :return: No return value.
        """

        rpkm_file = open(rpkm_file_path)
        matrix_file = open('../../../../Documents/Stanford/CS341_Data/donor_matrices_fixed/donor_' + donor + '.txt', 'w')

        firstLine = True
        donorColumns = donorsToColumns[donorId]
        for line in rpkm_file:
            if firstLine:
                firstLine = False
                continue
            else:
                lineAsArr = line.split('\t')
                for (j, rpkm_column) in enumerate(donorColumns):
                    expressionLevel = lineAsArr[rpkm_column]
                    if j == 0:
                        matrix_file.write(expressionLevel)
                    else:
                        matrix_file.write('\t' + expressionLevel)
            matrix_file.write('\n')

        matrix_file.close()
        rpkm_file.close()

    for (i, donor) in enumerate(donors):
        print 'Generating file for donor ', i, ': ', donor
        generate_donor_matrix(donor)
