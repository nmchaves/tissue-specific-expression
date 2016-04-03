# Convert RPKM Transcript Expression data to a matrix of
# tissue vs genes for a given user

# TODO: matrices are VERY sparse. Let's use scipy sparse matrices

#from scipy.sparse import lil_matrix
from scipy.sparse import csr_matrix
from scipy import io
import numpy as np


def saveStringArrayToFile(fName, array):
    """
        Saves the contents of |array| to the file specified by |fName|
    """
    f = open(fName, 'w')
    for elem in array[:-1]:
        f.write(elem + '\t')
    f.write(array[-1])

"""
TODO: open the GTEx_Data_V6_Annotations_SampleAttributesDS.txt file
and go through every sample (e.g. GTEX-111CU-1826-SM-5GZYN) to get its
tissue type. This field is called "SMTS" (see the GTEx_Data_V6_Annotations_SampleAttributesDD.xlsx
spreadsheet for descriptions of other fields). The SMTS field should be the 6th column
in the SampleAttributesDS.txt file, but to be safe I didn't hardcode this.

"""

attributesFile = open('../../../Downloads/GTEx_Data_V6_Annotations_SampleAttributesDS.txt')

tissueTypes = []
samplesToTissues = {}  # Dictionary mapping each sample ID to its tissue type
firstLine = True
for line in attributesFile:
    lineAsArr = line.split('\t')
    if firstLine:
        SAMPID_colno = lineAsArr.index('SAMPID')  # column of sample ID
        SMTS_colno = lineAsArr.index('SMTS')  # column of sample tissue type
        firstLine = False
    else:
        sampleId = lineAsArr[SAMPID_colno]
        tissueType = lineAsArr[SMTS_colno]
        if tissueType not in tissueTypes:
            tissueTypes.append(tissueType)
        samplesToTissues[sampleId] = tissueType

"""
The sample ID for an RNA-Seq or genotype sample is made up of the following 3 components separated by a dash, as exemplified with the example "GTEX-14753-1626-SM-5NQ9L":

1) "GTEX-YYYYY" (e.g. GTEX-14753) represents the GTEx donor ID. This ID should be used to link between the various RNA-Seq and genotype samples that come from the same donor.

2)"YYYY" (e.g., "1626") mostly refers to the tissue site, BUT we do not recommend using it for tissue site designation. Sometimes sample mix-ups occur,
and will be corrected however this part of the ID will not change when that happens. The accurate tissue site designation for all samples can be obtained from the
"Tissue Site Detail field" (encoded as "SMTSD") in the Sample Attributes file [Datasets->Download->GTEx_Data_V6_Annotations_SampleAttributesDS.txt].

3) "SM-YYYYY" (e.g., SM-5NQ9L) is the RNA or DNA aliquot ID used for sequencing.
'Y' stands for any number or capital letter.
"""

expressionFile = open('../../../Downloads/GTEx_Analysis_v6_RNA-seq_Flux1.6_transcript_rpkm.txt')

"""
 There are 195748 lines in the transcript_rpkm file

 There are a total of 8559 columns in the RPKM dataset. The 1st 4 columns are:
 1) TargetID
 2) Gene_Symbol
 3) Chr
 4) Coord

 The next 8555 columns are expression levels for a given user/tissue.
"""

columnsToSamples = []  # columnsToSamples[column] = id of the sample at this column in RPKM file
donorData = {}  # each element of this dictionary should be a list of tuples of the form {gene, tissue, exp level}
                # for a given donor
genes = []
firstLine = True
i = 0
for line in expressionFile:
    print i
    i += 1
    #if i > 1000:
    #    break
    lineAsArr = line.split("\t")
    if firstLine:
        for (col, sampleId) in enumerate(lineAsArr[4:-1]):
            columnsToSamples.append(sampleId)
            donorId = sampleId.split('-')[1]
            donorData[donorId] = []
        sampleId = (lineAsArr[-1])[:-1] # need to remove newline from the last sample ID
        columnsToSamples.append(sampleId)
        donorId = sampleId.split('-')[1]
        donorData[donorId] = []
        firstLine = False
    else:
        targetId = lineAsArr[0]
        geneSymbol = lineAsArr[1]
        if geneSymbol not in genes:
            genes.append(geneSymbol)
        chr = lineAsArr[2]
        coord = lineAsArr[3]  # pretty sure this is the genomic coordinate
        for (col, expLevel) in enumerate(lineAsArr[4:]):
            if expLevel > 0:
                sampleId = columnsToSamples[col]
                donorId = sampleId.split('-')[1]
                if donorId == '111CU':
                    tissueType = samplesToTissues[sampleId]
                    donorData[donorId].append((geneSymbol, tissueType, expLevel))

# Sort the genes and tissues, then generate matrices
genes = sorted(genes)
geneIndexes = {}
for (i, gene) in enumerate(genes):
    geneIndexes[gene] = i

tissueTypes = sorted(tissueTypes)
tissueIndexes = {}
for (i, tissue) in enumerate(tissueTypes):
    tissueIndexes[tissue] = i

#for donor in donorData.keys():
#    if donor == '111CU':
donor = '111CU'

m = np.zeros((len(genes), len(tissueTypes)))
#m = lil_matrix((len(genes), len(tissueTypes)))  # lil sparse format is good for
                                                        # appending to sparse matrices

fileName = 'genes_vs_tissues_2_' + donor
dataList = donorData[donor]
for dataTuple in dataList:
    row = geneIndexes[dataTuple[0]]  # Index of gene symbol in the sorted order
    col = tissueIndexes[dataTuple[1]]  # Index of tissue type in the sorted order
    m[row, col] = dataTuple[2]  # Expression Level

io.mmwrite(fileName, csr_matrix(m))  # Save in sparse matrix format
#np.savetxt(fileName, m, delimiter='\t')

'''
if tissueTypes[0] == '':
    # NOTE: there is some empty string tissue type, so the 1st element
    # in this array is an empty string
    tissueTypes[0] = 'BLANK'
saveStringArrayToFile('tissues.txt', tissueTypes)
saveStringArrayToFile('genes.txt', genes)
'''

#carray = np.asarray(tissueTypes)
#np.savetxt('tissues.txt', carray, delimiter='\t')
#print donorData['111CU']
"""
 Tissue type is SMTS
"""