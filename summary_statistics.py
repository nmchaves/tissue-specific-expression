"""
    This file summarizes the expression level of a GTEx RPKM file.

    If you want to run this script yourself, make sure to edit the file
    paths inside main.

    This script computes the following summary statistics

"""

# from scipy import stats
import numpy as np
import pandas as pd

"""
*********************
        Main
*********************
"""

def short_summary(filename):
    # compute and store the mean expression levels for each targetID
    meanExpLevels = [] 
    rpkm_file = open(filename)
    firstTarget = True
    firstLine = True
    for line in rpkm_file:
        if firstLine:
            firstLine = False
            continue
        else:
            lineAsArr = line.split('\t')
            expLevels = [float(i)  for i in lineAsArr]
            meanExpLevels.append(sum(expLevels))
    rpkm_file.close()

    # print meanExpLevels
    series = pd.Series(meanExpLevels)
    print series.describe()
    

if __name__ == "__main__":

    filename = 'tissue_test.txt' 
    short_summary(filename)

