"""
    This file performs the preprocessing of the GTEx RPKM file:
    "GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt"

    Running this script generates 2 files:

    1) donors.txt: This file will contain a tab-separated list of all
    donor ID's in the rpkm file (there will be 1 donor per column in this file)

    2) nonzeroTargetIds.txt: This file will contain each of the targetIds
    in the rpkm file that has at least 1 sample with a nonzero
    expression level.
"""

def generateDonorsFile(samples):
    """
    This function generates a tab-separated list of all
    donor ID's in the rpkm file.

    :param samples: An array containing all sample IDs in the rpkm file.
    :return: No return value.
    """

    donorFile = open('donors.txt', 'w')
    donorsSeen = {}
    firstDonor = True
    for sampleId in samples:
        donorId = sampleId.split('-')[1]
        if donorId not in donorsSeen:
            donorsSeen[donorId] = True
            if firstDonor:
                donorFile.write(donorId)
                firstDonor = False
            else:
                donorFile.write('\t' + donorId)


"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    path_to_rpkm_file = '../../../Downloads/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt'
    rpkm_file = open(path_to_rpkm_file)

    targetIds_file = open('nonzeroTargetIds.txt', 'w')
    targetIdsWritten = {}   # Dictionary containing all targetIds that have already been written
    firstTarget = True

    firstLine = True
    startIndexOfSamples = 4   # The 1st 4 columns of the 1st row of the RPKM file are not sample IDs.
    for line in rpkm_file:
        lineAsArr = line.split('\t')
        if firstLine:
            generateDonorsFile(lineAsArr[4:])
            firstLine = False
            continue
        else:
            # Generate the nonzero targetId file
            targetId = lineAsArr[0]
            if targetId not in targetIdsWritten:
                for expLevel in lineAsArr[startIndexOfSamples:]:
                    if float(expLevel) > 0:
                        if firstTarget:
                            targetIds_file.write(targetId)
                            firstTarget = False
                        else:
                            targetIds_file.write('\t' + targetId)
                        targetIdsWritten[targetId] = True
                        break

