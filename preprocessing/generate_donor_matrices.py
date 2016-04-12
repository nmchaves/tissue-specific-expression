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
        break # Only examine the 1st line of rpkm file

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

    donors = getArrayFromFile('donors.txt')
    targetIdsToSkip = getArrayFromFile('zeroTargetIds.txt')
    nonzeroTargetIds = getArrayFromFile('nonzeroTargetIds.txt')
    numTargetIds = len(targetIdsToSkip) + len(nonzeroTargetIds)

    rpkm_file_path = '../../../Downloads/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt'
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
        matrix_file = open('donor_matrices/donor_' + donor + '.txt', 'w')

        firstLine = True
        row = 0
        donorColumns = donorsToColumns[donorId]
        for line in rpkm_file:
            #if row > 1000:
            #    break  # TODO: remove this! it's just for testing
            if firstLine:
                firstLine = False
                continue
            else:
                targetId = line[0:line.index('\t')]

                if targetId in targetIdsToSkip:
                    # Write a row full of zeros
                    matrix_file.write('0')
                    for i in range(len(donorColumns)-1):
                        matrix_file.write('\t0')
                else:
                    lineAsArr = line.split('\t')
                    for (j, rpkm_column) in enumerate(donorColumns):
                        expressionLevel = lineAsArr[rpkm_column]
                        if j==0:
                            matrix_file.write(expressionLevel)
                        else:
                            matrix_file.write('\t' + expressionLevel)
                matrix_file.write('\n')
                row += 1

        matrix_file.close()
        rpkm_file.close()

    for (i, donor) in enumerate(donors):
        print 'donor ', i, ': ', donor
        generate_donor_matrix(donor)
        break # TODO: remove this break. it's just for testing
