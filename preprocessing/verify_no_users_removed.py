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

    NUM_TRANSCRIPTS = 10000  # Adjust this if not using 10,000 transcripts

    # If you run this yourself, you should set these paths properly!!
    donors_file_path = '../data/donors.txt'
    donor_data_path = '../../../../Documents/Stanford/CS341_Data/donor_matrices_fixed/'
    donor_file_prefix = 'donor_'
    donors = getArrayFromFile(donors_file_path)

    #output_file = open('../data/donors_mean_exp_levels_after_filtering.txt', 'w')
    #output_file.write('Donor_ID\tNum_Total_Samples\tNum_Unusable_Samples\tOverall_Mean_Exp_Level\tMean_Exp_Level_List\n')

    mean_expression_levels = []  # Array containing mean expression level for each user
    print 'numDonors: ', len(donors)

    for donor in donors:
        if donor == 'ZXG5':
            continue
        file_name = donor_data_path + donor_file_prefix + donor + '.txt'
        donor_file = open(file_name)
        levels_so_far = []
        firstRow = True
        num = 0
        for line in donor_file:
            num += 1
            exp_levels = line.split('\t')
            if firstRow:
                for exp_level in exp_levels:
                    levels_so_far.append(float(exp_level))
                firstRow = False
            else:
                for (i, exp_level) in enumerate(exp_levels):
                    try:
                        levels_so_far[i] += float(exp_level)
                    except:
                        print donor, ' line: ', num
                        print 'issue: ', exp_level

        num_samples = len(levels_so_far)
        mean_exp_by_sample = [1.0 * level / NUM_TRANSCRIPTS for level in levels_so_far]
        for (i, mean_exp) in enumerate(mean_exp_by_sample):
            if mean_exp < 0.01:
                print 'exp_level = ', mean_exp, ' for donor: ', donor, ' in column: ', i
        mean_over_all_samples = sum(mean_exp_by_sample)/len(mean_exp_by_sample)
        mean_expression_levels.append(mean_over_all_samples)
        donor_file.close()

        num_unusable_samples = 0
        #output_file.write(donor + '\t' + str(num_samples) + '\t' + str(num_unusable_samples) + '\t'
        #                  + str(mean_over_all_tissues) + '\t' + '\n')

    #output_file.close()
    # Sort by expression level
