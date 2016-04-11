"""
    This script generates a file containing a tab-delimited matrix
    for each tissue.

    The matrix files are generated incrementally. If this program gets interrupted at some
    point, then we'll have saved the files for each tissue that has been processed so far.

"""

def getTissueColumns(tissue, path_to_rpkm_file):
    """
    Get the columns of the RPKM file with a sample from this tissue.

    :param path_to_rpkm_file: Path to the RPKM file
    :return: LIst of columns
    """
    columns = []
    # Get a dictionary whose keys are sample IDs of this tissue's samples.
    samples = getDictFromFile('../../../Documents/Stanford/CS341_Data/tissue_metadata/tissue_meta_' + tissue + '.txt')
    rpkm_file = open(path_to_rpkm_file)
    for line in rpkm_file:
        lineAsArr = line.split('\t')
        for (col, sampleId) in enumerate(lineAsArr):
            if col >= 4:
                if sampleId in samples:
                    columns.append(col)
        break  # Only examine the 1st line of rpkm file
    return columns

def getDictFromFile(path):
    """
    This function takes a path to a file containing a 1-line tab-delimited
    text file and converts that file into a dictionary.

    :param path: Path to the file
    :return: The text file
    """
    f = open(path)
    dictObject = {}
    for line in f:
        lineAsArr = line.split('\t')
        break
    for element in lineAsArr:
        dictObject[element] = True
    return dictObject


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

    tissues = getArrayFromFile('tissues.txt')
    targetIdsToSkip = getArrayFromFile('zeroTargetIds.txt')
    nonzeroTargetIds = getArrayFromFile('nonzeroTargetIds.txt')
    numTargetIds = len(targetIdsToSkip) + len(nonzeroTargetIds)

    rpkm_file_path = '../../../Documents/Stanford/CS341_Data/transcript_rpkm_top_10000_var.txt'

    def generate_tissue_matrix(tissue):
        """
        This function generates the tissue matrix for |tissue|.
        After this function completes there will be a file containing
        a tab-delimited matrix. The matrix columns correspond to samples and the
        matrix rows correspond to target IDs. The cell values are expression
        levels.

        :param tissue: Name of the tissue
        :return: No return value.
        """

        rpkm_file = open(rpkm_file_path)
        matrix_file = open('tissue_matrices/tissue_' + tissue + '.txt', 'w')

        firstLine = True
        row = 0
        tissueColumns = getTissueColumns(tissue, rpkm_file_path)
        for line in rpkm_file:
            if firstLine:
                firstLine = False
                continue
            else:
                targetId = line[0:line.index('\t')]

                if targetId in targetIdsToSkip:
                    # Write a row full of zeros
                    matrix_file.write('0')
                    for i in range(len(tissueColumns)-1):
                        matrix_file.write('\t0')
                else:
                    lineAsArr = line.split('\t')
                    for (j, rpkm_column) in enumerate(tissueColumns):
                        expressionLevel = lineAsArr[rpkm_column]
                        if j == 0:
                            matrix_file.write(expressionLevel)
                        else:
                            matrix_file.write('\t' + expressionLevel)
                matrix_file.write('\n')
                row += 1

        matrix_file.close()
        rpkm_file.close()

    for (i, tissue) in enumerate(tissues):
        print 'Generating file for tissue ', i, ': ', tissue
        generate_tissue_matrix(tissue)
