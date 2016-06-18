"""
    This file performs all of the preprocessing steps needed for the GTEx RPKM file:
    "GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt"

    If you want to run this script yourself, make sure to update the file
    paths inside main to reflect the paths on your own system.

    Running this script generates the following 8 types of files.

    1) donors.txt: Contains a tab-separated list of all
    donor ID's in the rpkm file (See generateDonorsFile())

    2) donorTissues.txt: Each row contains a tab-separted list of the form:
    donorID     number_of_tissues_from_this_donor   tissue1_name    tissue2_name...
    (See generateDonorTissuesFile())

    3) tissues.txt: Tab-delimited file of all SMTSD tissues. These are NOT all in alphabetical
        order. Rather, they are in order of their appearance in the RPKM file.

    4) donor_meta_[donor Id].txt: The meta file for each donor.
    Output file formatted as:

    Donor_ID    Gender     Age     DTHHRDY      Mean_Expression_Over_All_Tissues    Num_Samples
    ABCDEF      M          20      1            110.3                               10
    ----------
    Sample_ID   Tissue          Mean_Expression     Variance    25th percentile....
    ABCDEFG...  Brain - ...     98.76
    ABCDEFG...  Fat - ...       10.2.6

    5) tissue_meta_[tissue name].txt: Meta file for each tissue. This is a tab-delimited
    list of all sample IDs corresponding to this tissue

"""


import numpy as np

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
    rpkm_file.close()


def samplesToTissuesDict(attributes_path):
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
            if sampleId[-1] == '\n':
                sampleId = sampleId[0:-1]
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

    samplesToTissues = samplesToTissuesDict(attributes_path)

    donorsToTissues = {}  # maps a donor ID to list of tissues sampled
    numSamples = len(samples)
    for (i, sampleId) in enumerate(samples):
        if i == numSamples-1:
            sampleId = samples[-1][:-1] # need to remove newline from the last sample ID
        donorId = sampleId.split('-')[1]
        if donorId == 'OHPJ':
            print donorId
        tissue = samplesToTissues[sampleId]
        if donorId not in donorsToTissues:
            donorsToTissues[donorId] = [tissue]
        else:
            if tissue not in donorsToTissues[donorId]:
                donorsToTissues[donorId].append(tissue)

    # Go through dictionary and save the file
    donors_tissue_file = open('donorTissues.txt', 'w')
    for donorId in sorted(donorsToTissues.keys()):
        tissueList = donorsToTissues[donorId]
        numTissuesStr = str(len(tissueList))
        donors_tissue_file.write(donorId + '\t' + numTissuesStr)
        for tiss in tissueList:
            donors_tissue_file.write('\t' + tiss)
        donors_tissue_file.write('\n')
    donors_tissue_file.close()
    rpkm_file.close()


def generateTissuesFile(rpkm_path, attributes_path):

    samplesToTissues = samplesToTissuesDict(attributes_path)

    tissues = []
    rpkm_file = open(rpkm_path)
    for line in rpkm_file:
        samples = line.split('\t')[4:]
        break
    samples[-1] = samples[-1][:-1]  # Remove newline from last sample ID

    for sample in samples:
        tissue = samplesToTissues[sample]
        if tissue not in tissues:
            tissues.append(tissue)

    tissue_file = open('tissues.txt', 'w')
    tissue_file.write(tissues[0])
    for t in tissues[1:]:
        tissue_file.write('\t' + t)
    rpkm_file.close()


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
    rpkm_file.close()

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
    rpkm_file.close()

def buildDonorsToColumnsAndSamplesDict(path_to_rpkm_file):
    """
    This function creates a dictionary that maps a donor to the list
    of columns/sample IDs that correspond to this donor in the RPKM file.

    :param path_to_rpkm_file: Path to the RPKM txt file
    :return: The dictionary
    """

    rpkm_file = open(path_to_rpkm_file)

    donorsDict = {}
    for line in rpkm_file:
        lineAsArr = line.split('\t')
        for (col, sampleId) in enumerate(lineAsArr):
            if col >= 4:
                if sampleId[-1] == '\n':
                    sampleId = sampleId[0:-1]
                donorId = sampleId.split('-')[1]
                if donorId not in donorsDict:
                    donorsDict[donorId] = [(col, sampleId)]
                else:
                    donorsDict[donorId].append((col, sampleId))
        break  # Only examine the 1st line of rpkm file

    rpkm_file.close()
    return donorsDict

def buildDonorsToPhenotypeDict(path_to_phenotype_file):
    phen_file = open(path_to_phenotype_file)
    firstLine = True
    donorsToPhenotypeDict = {}
    for line in phen_file:
        if firstLine:
            firstLine = False
            continue
        else:
            vals = line.split('\t')
            donor = vals[0].split('-')[1]
            # Add the phenotypes. Must remove the newline from last phenotype value!
            donorsToPhenotypeDict[donor] = (vals[1], vals[2], vals[3][0:-1])
    phen_file.close()
    return donorsToPhenotypeDict


def generateDonorMetaFiles(rpkm_path, attributes_path):

    def getStatistics(values):
        mean = np.mean(values)
        min = np.min(values)
        std_dev = np.std(values)
        pct_25 = np.percentile(values, 25)
        pct_75 = np.percentile(values, 75)
        max = np.max(values)
        return str(mean) + '\t' + str(std_dev) + '\t' + str(min) + '\t' + \
               str(pct_25) + '\t' + str(pct_75) + '\t' + str(max)

    samplesToTissues = samplesToTissuesDict(attributes_path)
    donorsToColumnsAndSamples = buildDonorsToColumnsAndSamplesDict(rpkm_path)
    donorsToPhenotypes = buildDonorsToPhenotypeDict('../data/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt')

    donors = sorted(donorsToColumnsAndSamples.keys())
    for donor in donors:
        print donor
        columnsAndSamples = donorsToColumnsAndSamples[donor]

        sample_expressions = {}
        for _, sample in columnsAndSamples:
            sample_expressions[sample] = np.zeros(NUM_TRANSCRIPTS_TO_RETAIN)

        donor_mean_exp_level = 0
        rpkm_file = open(rpkm_path)
        firstLine = True
        i = 0
        for line in rpkm_file:
            if firstLine:
                firstLine = False
                continue
            vals = line.split('\t')
            for col, sample in columnsAndSamples:
                expression_level = float(vals[col])
                donor_mean_exp_level += expression_level
                sample_expressions[sample][i] = expression_level
            i += 1

        num_samples = len(columnsAndSamples)
        donor_mean_exp_level /= (NUM_TRANSCRIPTS_TO_RETAIN*num_samples)

        # Write out results to meta_file
        donor_meta_file = open('../data/Donor_Metadata_Enhanced/donor_meta_' + donor + '.txt', 'w')
        header = 'Donor_ID\tGender\tAge\tDTHHRDY\tMean_Expression_Over_All_Tissues\tNum_Samples\n'
        donor_meta_file.write(header)
        phenotype = donorsToPhenotypes[donor]
        gender = phenotype[0]
        age = phenotype[1]
        death_type = phenotype[2]
        donor_meta_file.write(donor + '\t' + gender + '\t' + age + '\t' +
                              death_type + '\t' + str(donor_mean_exp_level) + '\t' + str(num_samples) + '\n')
        donor_meta_file.write('----------\n')
        header2 = 'Sample_ID\tTissue\tMean_Expression\tStd_Dev\tMin\t25th_Pct\t75th_Pct\tMax\n'
        donor_meta_file.write(header2)
        for _, sample in columnsAndSamples:
            donor_meta_file.write(sample + '\t' + samplesToTissues[sample] + '\t')
            expression_levels = sample_expressions[sample]
            donor_meta_file.write(getStatistics(expression_levels) + '\n')
        donor_meta_file.close()
        rpkm_file.close()


def generateTissueMetaFiles(rpkm_path, attributes_path):

    tissues = getArrayFromFile('../data/tissues.txt')
    sampleToTissues = samplesToTissuesDict(attributes_path)
    tissueToSamples = {}

    for tissue in tissues:
        tissueToSamples[tissue] = []
        tissue_metafile = open('../data/tissue_metadata/tissue_meta_' + tissue + '.txt', 'w')
        tissue_metafile.write('# Sample Index\tSample ID\n')
        rpkm_file = open(rpkm_path)
        for line in rpkm_file:
            line = line.rstrip()
            sampleIds = line.split('\t')[4:]
            break
        #sampleIds[-1] = sampleIds[-1][:-1]

        for (i, sampleId) in enumerate(sampleIds):
            sampleTissue = sampleToTissues[sampleId]
            if sampleTissue == tissue:
                tissue_metafile.write(str(i) + '\t' + sampleId + '\n')

        tissue_metafile.close()
        rpkm_file.close()



def filter_by_go_annotation(path_to_rpkm_file, path_to_output_file, path_to_list):
    """
    This function goes through each line of the rpkm file. If the current transcript
    corresponds to a gene in Gene Ontology, then the current line is written out
    to the file "transcript_rpkm_in_go.txt".

    :param path_to_rpkm_file: input file name
    :param path_to_output_file: ouput file name
    :param path_to_list: gene conversion file name
    :return: No return value.
    """

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

    firstLine = True
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            # Write the header row
            rpkm_file_retained.write(line)
            continue
        else:
            transcriptId = line[0:line.index('\t')]
            if '.' in transcriptId:
                transcriptId = transcriptId[0:transcriptId.index('.')]

            if transcriptId in transcriptIdsInGO:
                rpkm_file_retained.write(line)

    rpkm_file_retained.close()
    return


def filter_remove_zero_expression(path_to_rpkm_file, path_to_output_file):
    rpkm_file = open(path_to_rpkm_file)
    new_rpkm_file = open(path_to_output_file, 'w')

    firstLine = True
    num_rows_retained = 0
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            new_rpkm_file.write(line)
            continue

        exp_levels = line.rstrip().split('\t')[4:]
        nonzero = False
        for exp_level in exp_levels:
            if float(exp_level) != 0:
                nonzero = True
                break

        if nonzero:
            num_rows_retained += 1
            new_rpkm_file.write(line)

    rpkm_file.close()
    new_rpkm_file.close()
    print "Number of transcripts retained: ", num_rows_retained


def filter_by_var(path_to_rpkm_file, path_to_output_file, NUM_TRANSCRIPTS_TO_RETAIN):
    """
        This function should be ran AFTER filter_by_go.py

        This script takes the top |NUM_TRANSCRIPTS_TO_RETAIN| rows by variance and
        saves them into a text file. In other words, this file removes the low variance
        rows from the rpkm file's matrix.

        :param path_to_rpkm_file: input file name
        :param path_to_output_file: ouput file name
        :param NUM_TRANSCRIPTS_TO_RETAIN: number of transcripts to retain
        :return: No return value.

    """

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
    return



def getArrayFromFile(path):
    """
    This function takes a path to a file containing a tab-delimited
    text file and converts that file into an array.

    :param path: Path to the file
    :return: The array
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

    raw_rpkm_file = '../../../../Downloads/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt'
    path_to_attributes_file = '../../../../Downloads/GTEx_Data_V6_Annotations_SampleAttributesDS.txt'
    NUM_TRANSCRIPTS_TO_RETAIN = 10000  # For certain purposes (eg t-SNE, only retain 10,000 transcripts

    # Operate on the RPKM expression levels matrix itself.
    local_data_dir = '../../CS341_Data/'  # Change this to reflect your local machine.
    data_dir = '../data/'
    geneListGO = data_dir + 'trans_gene_name_filtered_by_GO.txt'
    output_go  = local_data_dir + 'transcript_rpkm_in_go.txt'
    output_go_nonzero = local_data_dir + 'transcript_rpkm_in_go_nonzero.txt'
    output_var = local_data_dir + 'transcript_rpkm_in_go_top_' + str(NUM_TRANSCRIPTS_TO_RETAIN) + '_var.txt'

    print 'Filter 1: Filtering by GO gene list'
    filter_by_go_annotation(raw_input, output_go, geneListGO)
    print 'Filter 2: Removing transcripts with all 0 expression'
    filter_remove_zero_expression(output_go, output_go_nonzero)
    print 'Filter 3: Filtering by variance'
    filter_by_var(output_go_nonzero, output_var, NUM_TRANSCRIPTS_TO_RETAIN)

    # Generate MetaData Files
    generateDonorsFile(raw_rpkm_file)
    generateDonorMetaFiles(output_var, path_to_attributes_file) # TODO: run this on the full RPKM matrix
    generateDonorTissuesFile(raw_rpkm_file, path_to_attributes_file)
    generateTissueMetaFiles(raw_rpkm_file, path_to_attributes_file)



