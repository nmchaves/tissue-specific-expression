"""
    This file performs the preprocessing of the GTEx RPKM file:
    "GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt"

    If you want to run this script yourself, make sure to edit the file
    paths inside main.

    Running this script generates the following 5 files (the 5th "file" is actually
    many files--1 for each donor):

    1) donors.txt: This file will contain a tab-separated list of all
    donor ID's in the rpkm file (there will be 1 donor per column in this file)
    (See generateDonorsFile())

    2) donorTissues.txt: Each row contains a tab-separted list of the form:
    donorID     number_of_tissues_from_this_donor   tissue1_name    tissue2_name...
    (See  generateDonorTissuesFile())

    3) nonzeroTargetIds.txt: This file will contain each of the targetIds
    in the rpkm file that has at least 1 sample with a nonzero
    expression level.
    (See generateTargetIdFiles())

    4) zeroTargetIds.txt: This file will contain each of the targetIds
    in the rpkm file that has all 0 expression levels.
    (See generateTargetIdFiles())

    5) donor_meta_[donor Id].txt: The meta file for each donor. Column 1 is sample
    IDs. Column 2 is tissue type

"""

def generateDonorsFile(path):
    """
    This function generates a file containing a tab-separated
    list of all donor ID's in the rpkm file.

    :param path: Path to the rpkm file.
    :return: No return value.
    """
    rpkm_file = open(path)
    for line in rpkm_file:
        samples = line.split('\t')[4:]
        break

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
    donorFile.close()


def generateSamplesToTissuesDict(attributes_path):
    attributes_file = open(attributes_path)
    tissueTypes = []
    samplesToTissues = {}  # Dictionary mapping each sample ID to its tissue type
    firstLine = True
    for line in attributes_file:
        lineAsArr = line.split('\t')
        if firstLine:
            SAMPID_colno = lineAsArr.index('SAMPID')  # column of sample ID
            SMTS_colno = lineAsArr.index('SMTSD')  # column of sample tissue type
            firstLine = False
        else:
            if lineAsArr[0][0:4] != 'GTEX':
                continue
            sampleId = lineAsArr[SAMPID_colno]
            tissueType = lineAsArr[SMTS_colno]
            if tissueType not in tissueTypes:
                tissueTypes.append(tissueType)
            samplesToTissues[sampleId] = tissueType
    return samplesToTissues

def generateDonorTissuesFile(rpkm_path, attributes_path):

    rpkm_file = open(rpkm_path)

    for line in rpkm_file:
        samples = line.split('\t')[4:]
        break

    samplesToTissuesDict = generateSamplesToTissuesDict(attributes_path)
    donorsToTissues = {}  # maps a donor ID to list of tissues sampled
    numSamples = len(samples)
    for (i, sampleId) in enumerate(samples):
        if i == numSamples-1:
            sampleId = samples[-1][:-1] # need to remove newline from the last sample ID
        donorId = sampleId.split('-')[1]
        tissue = samplesToTissuesDict[sampleId]
        if donorId not in donorsToTissues:
            donorsToTissues[donorId] = [tissue]
        else:
            if tissue not in donorsToTissues[donorId]:
                donorsToTissues[donorId].append(tissue)

    # Go through dictionary and save the file
    donors_tissue_file = open('donorTissues.txt', 'w')
    for donorId in sorted(donorsToTissues):
        tissueList = donorsToTissues[donorId]
        numTissuesStr = str(len(tissueList))
        donors_tissue_file.write(donorId + '\t' + numTissuesStr)
        for tiss in tissueList:
            donors_tissue_file.write('\t' + tiss)
        donors_tissue_file.write('\n')
    donors_tissue_file.close()



def generateTargetIdFiles(path):
    rpkm_file = open(path)

    targetIds_file = open('nonzeroTargetIds.txt', 'w')
    targetIdsWritten = {}   # Dictionary containing all targetIds that have already been written
                            # out to nonzeroTargetIds.txt


    # Generate targetIds file
    firstTarget = True
    firstLine = True
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            continue
        else:
            targetId = line[0:line.index('\t')]
            if targetId not in targetIdsWritten:
                lineAsArr = line.split('\t')
                for expLevel in lineAsArr[4:]:
                    if float(expLevel) > 0:
                        if firstTarget:
                            targetIds_file.write(targetId)
                            firstTarget = False
                        else:
                            targetIds_file.write('\t' + targetId)
                        targetIdsWritten[targetId] = True
                        break  # go to next line in file
    targetIds_file.close()

    # Generate zero targetIds file
    zeroTargetIds_file = open('zeroTargetIds.txt', 'w')
    rpkm_file = open(path)
    firstLine = True
    firstTarget = True
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            continue
        else:
            targetId = line[0:line.index('\t')]
            if targetId not in targetIdsWritten:
                if firstTarget:
                    zeroTargetIds_file.write(targetId)
                    firstTarget = False
                else:
                    zeroTargetIds_file.write('\t' + targetId)
    zeroTargetIds_file.close()

def generateDonorMetaFiles(rpkm_path, attributes_path):

    samplesToTissues = generateSamplesToTissuesDict(attributes_path)

    rpkm_file = open(rpkm_path)
    for line in rpkm_file:
        samples = line.split('\t')[4:]
        break
    samples[-1] = samples[-1][:-1]  # Remove the newline character from the last sample ID
    rpkm_file.close()

    donorMetaData = {}
    for sampleId in samples:
        donorId = sampleId.split('-')[1]
        tissue = samplesToTissues[sampleId]
        if donorId not in donorMetaData:
            donorMetaData[donorId] = [(sampleId, tissue)]
        else:
            donorMetaData[donorId].append((sampleId, tissue))

    for donorId in donorMetaData.keys():
        donor_meta_file = open('donor_matrices/donor_meta_' + donorId + '.txt', 'w')
        metaData = donorMetaData[donorId]
        for (sampleId, tissue) in metaData:
            donor_meta_file.write(sampleId + '\t' + tissue + '\n')
        donor_meta_file.close()


"""
*********************
        Main
*********************
"""
if __name__ == "__main__":

    path_to_rpkm_file = '../../../Downloads/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt'
    path_to_attributes_file = '../../../Downloads/GTEx_Data_V6_Annotations_SampleAttributesDS.txt'

    generateDonorsFile(path_to_rpkm_file)
    generateDonorMetaFiles(path_to_rpkm_file, path_to_attributes_file)
    generateDonorTissuesFile(path_to_rpkm_file, path_to_attributes_file)
    generateTargetIdFiles(path_to_rpkm_file)


