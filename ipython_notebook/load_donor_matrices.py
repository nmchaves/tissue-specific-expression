
# coding: utf-8

# In[9]:

# get_ipython().magic(u'matplotlib inline')
# from sklearn import decomposition
# from sklearn import manifold
# import matplotlib.pyplot as plt


# In[ ]:

def construct_matrix(spec_donors=[],reject_donors=[],donor_sex=None,donor_age=[],spec_tissues=[],reject_tissues=[],
                       top_tissues=None,min_tissue_count=None):
    #     This function builds a numpy matrix and associated data arrays for the specified
    
    #     --Parameters--
    #     spec_donors: Array of donor ID strings to include. Leave blank for all donors.
    #     reject_donors: Array of donor IDs to remove. Leave blank to reject no donors.
    #     donor_sex: '1' for only males, '2' for only females, None (or blank) for all donors.
    #     donor_age: Array of age-decade strings to include, e.g. ['20','30','40']. Leave blank for all donors.
    #     spec_tissues: Array of tissue names to include. Leave blank for all tissues.
    #     reject_tissues: Array of tissue names to remove. Leave blank to reject no tissues.
    #     top_tissues: Number of most commonly sampled tissues to include. Leave blank to include all tissues.
    #     min_tissue_count: Minimum number of samples to include a tissue. Leave blank to include all tissues.

    #     Note: spec_donors and reject_donors are mutually exclusive.
    #     Note: spec_tissues and reject_tissues are mutually exclusive.
    #     Note: top_tissues and min_tissue_count are mutually exclusive.
    
    #     --Returns--
    #     multi_matrix: Numpy matrix of 10,000 rows, column for each included sample
    #     point_ID: Array of donor IDs for each column in multi_matrix
    #     point_sex: Array of donor sex for each column in multi_matrix
    #     point_age: Array of donor age for each column in multi_matrix
    #     point_tissue: Array of tissue type for each column in multi_matrix
    #     point_tissue_group: Array of tissue group (just first term of tissue type) for each column in multi_matrix

    # DONORS BY ID
    if not spec_donors:
        all_donors = open('../data/donors_list.txt')
        donor_list = [ID[0:-1] for ID in all_donors if ID[0:-1] not in reject_donors]
        all_donors.close()
    else:
        donor_list = spec_donors
    # dictionary of donor IDs, and an array that will be [sex,age]
    donor_dict = dict((ID,[]) for ID in donor_list)

    # DONORS BY AGE AND SEX
    donor_info = open('../data/donor_info.txt')
    for line in donor_info:
        # info is [ID,sex,age,death]
        info = line.split('\t')
        if info[0] in donor_list:
            # check sex
            if donor_sex and info[1] != donor_sex:
                del donor_dict[info[0]]
                continue
            else:
                donor_dict[info[0]].append(info[1])
            # check age    
            age = info[2].split('-')[0]
            if donor_age and age not in donor_age:
                del donor_dict[info[0]]
                continue
            else:
                donor_dict[info[0]].append(age)
    donor_info.close()

    # TISSUES BY TYPE AND SAMPLE COUNT
    tissues = Counter()
    for ID in donor_dict.keys():
        metafile = open('../data/Donor_Metadata_Enhanced/donor_meta_'+ID+'.txt')
        headerLines = 4
        lineCounter = 0
        for line in metafile:
            # skip the first four lines of header
            if lineCounter < headerLines:
                lineCounter += 1
                continue 
            # look for tissue type listed in meta file
            tissue = line.split('\t')[1][0:-1]
            if spec_tissues and tissue in spec_tissues:
                tissues[tissue] = tissues.get(tissue,0) + 1
            elif not spec_tissues and tissue not in reject_tissues:
                tissues[tissue] = tissues.get(tissue,0) + 1
        metafile.close()
    if min_tissue_count:
        tissue_list = [key for key,value in tissues.iteritems() if value >= min_tissue_count]
    else:
        tissue_list = [key for key,value in tissues.most_common(top_tissues)]

    # CONSTRUCT MATRIX
    # initialize column (to be removed) and info types
    multi_matrix = np.zeros((10000,1))
    point_ID = []
    point_sex = []
    point_age = []
    point_tissue = []
    point_tissue_group = []
    # metadata on relevant points
    for ID in donor_dict.keys():
        metafile = open('../data/Donor_Metadata_Enhanced/donor_meta_'+ID+'.txt')
        # column indices for relevant tissues
        columns = []
        column = 0
        headerLines = 4
        lineCounter = 0
        for line in metafile:
            # deal with the first four lines of header
            if lineCounter < headerLines:
                lineCounter += 1
                continue 
            tissue = line.split('\t')[1][0:-1]      
            if tissue in tissue_list:
                columns.append(column)
                point_ID.append(ID)
                point_sex.append(donor_dict[ID][0])
                point_age.append(donor_dict[ID][1])
                point_tissue.append(tissue)
                point_tissue_group.append(tissue.split('-')[0])
            column = column + 1
        metafile.close()
        # get data
        donor_matrix = np.zeros((10000,len(columns)))
        row = 0
        donorfile = open('../data/donor_matrices_fixed/donor_'+ID+'.txt')
        for line in donorfile:
            # from file, take desired tissue columns and put in donor matrix
            values = [line.split('\t')[ind] for ind in columns]
            donor_matrix[row,:] = values
            row = row+1
        donorfile.close()
        # concatenate donor matrices
        multi_matrix = np.concatenate((multi_matrix,donor_matrix),axis=1)
    multi_matrix = np.delete(multi_matrix,0,1)
    
    print 'Matrix constructed with ' + str(multi_matrix.shape[1]) + ' samples!'
    return [multi_matrix, point_ID, point_sex, point_age, point_tissue, point_tissue_group]


# In[ ]:

# multi_matrix, point_ID, point_sex, point_age, point_tissue, point_tissue_group = construct_matrix()


# In[ ]:



