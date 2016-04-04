"""
    This script generates a file containing a sparse matrix (mtx extension)
    for each donor.

    The matrix files are generated incrementally. If this program is interrupted at some
    point, then we'll have saved the files for each donor that has been processed so far.

"""

from scipy.sparse import csr_matrix
from scipy import io
import numpy as np

def buildDonorsToColumnsDict(path_to_rpkm_file):
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

        break

    return donorsDict


def getArrayFromFile(path):
    f = open(path)
    for line in f:
        return line.split('\t')

def getDonorTissues():
    donor_tissues_file = open('donorTissues.txt')
    donorToTissueCountDict = {}
    for line in donor_tissues_file:
        lineAsArr = line.split('\t')
        donorToTissueCountDict[lineAsArr[0]] = int(lineAsArr[1])
    donor_tissues_file.close()
    return donorToTissueCountDict

"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    #tissues = getArrayFromFile('tissues.txt')
    donors = getArrayFromFile('donors.txt')
    targetIdsToSkip = getArrayFromFile('zeroTargetIds.txt')
    nonzeroTargetIds = getArrayFromFile('nonzeroTargetIds.txt')
    numTargetIds = len(targetIdsToSkip) + len(nonzeroTargetIds)

    rpkm_file_path = '../../../Downloads/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt'
    donorsToColumns = buildDonorsToColumnsDict(rpkm_file_path)

    donorToTissueCount = getDonorToTissueCounts()

    def generate_donor_matrix(donorId):
        '''
        This function generates the donor matrix for |donorId|.
        After this function completes there will be a file containing
        a sparse matrix. The matrix rows correspond to samples and the
        matrix columns correspond to target IDs. The cell values are expression
        levels.

        :param donorId: The donor's ID
        :return: No return value.
        '''

        numTissuesCurDonor = donorToTissueCount[donorId]
        print numTissuesCurDonor, ' ', numTargetIds
        m = np.zeros(shape=(numTargetIds, numTissuesCurDonor))
        rpkm_file = open(rpkm_file_path)

        firstLine = True
        row = 0
        for line in rpkm_file:
            if row > 1000:
                break  # TODO: remove this! it's just for testing
            if firstLine:
                firstLine = False
                continue
            else:
                targetId = line[0:line.index('\t')]

                if targetId in targetIdsToSkip:
                    continue
                else:
                    lineAsArr = line.split('\t')
                    for (j, rpkm_column) in enumerate(donorsToColumns[donorId]):
                        expLevel = lineAsArr[rpkm_column]
                        m[row, j] = expLevel
                row += 1

        file_name = 'donor_matrices/donor_' + donor
        io.mmwrite(file_name, csr_matrix(m))
        #cur_donor_metafile = open('donor_matrices/donor_meta_' + donor + '.txt', 'w')
        rpkm_file.close()

    for (i, donor) in enumerate(donors):
        if i==0:
            continue
        print i , 'th donor: ', donor
        generate_donor_matrix(donor)
        break # TODO: remove this break. it's just for testing
