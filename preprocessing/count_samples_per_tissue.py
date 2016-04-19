

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

    output_file = open('../data/tissue_counts_male_female.txt', 'w')
    output_file.write('Tissue_Name\tTotal_Samples\tMale_Samples\tFemale_Samples\n')
    tissues = getArrayFromFile('../data/tissues.txt')
    tissue_metadata_path = '../../../../Documents/Stanford/CS341_Data/tissue_metadata/tissue_meta_'

    donors_phenotype_file = open('../../../../Downloads/GTEx_Data_V6_Annotations_SubjectPhenotypesDS.txt')
    id_to_gender = {}
    # Convention in GTEX phenotypes file is male=1, female=2
    MALE = 1
    FEMALE = 2
    firstLine = True
    for line in donors_phenotype_file:
        if firstLine:
            firstLine = False
            continue
        else:
            lineAsArr = line.split('\t')
            donor_id = (lineAsArr[0].split('-'))[1]
            id_to_gender[donor_id] = int(lineAsArr[1])
    donors_phenotype_file.close()

    for tissue in tissues:
        output_file.write(tissue + '\t')
        tissue_meta_file = open(tissue_metadata_path + tissue + '.txt')
        tissue_sample_ids = []
        for line in tissue_meta_file:
            tissue_sample_ids = line.split('\t')
            break
        num_samples = len(tissue_sample_ids)
        output_file.write(str(num_samples) + '\t')

        num_males = 0
        num_females = 0
        for sample in tissue_sample_ids:
            donor_id = sample.split('-')[1]
            if id_to_gender[donor_id] == MALE:
                num_males += 1
            elif id_to_gender[donor_id] == FEMALE:
                num_females += 1
            else:
                print 'error!'
        output_file.write(str(num_males) + '\t')
        output_file.write(str(num_females) + '\n')
        tissue_meta_file.close()

    output_file.close()
