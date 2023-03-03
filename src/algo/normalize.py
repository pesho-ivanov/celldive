# A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis (2013)


import pandas as pd
import numpy as np
#import rpy2.robjects as robj
#from rpy2.robjects import numpy2ri
#from rpy2.robjects import pandas2ri
#from rpy2.robjects.packages import importr

import matplotlib.pyplot as plt
from utils import *

#def NormalizeDESeq(countTable):
#    deseq = importr("DESeq")
#
#    numpy2ri.activate()
#    pandas2ri.activate()
#    
#    # genes X cells
#    
#    #display(countTable.loc['290'])
#    #display(countTable['ENST00000390455'])
#    
#    #countTable['ENST00000390393'] = 100 
#    conditions = robj.FactorVector( 'c' * len(countTable.index) )
#    cds = deseq.newCountDataSet( countTable.T, conditions )
#    data = cds.slots['assayData']['counts']
#    #print(data)
#    sizeFactors = deseq.estimateSizeFactorsForMatrix( data )
#    #print('normalization sizeFactors ({}): {}'.format(len(list(sizeFactors)), list(sizeFactors)))
#    #return countTable / sizeFactors
#    return countTable.div(sizeFactors, axis='index')

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def NormalizeByTotalSum(countTable):
    #print('NormalizeByTotalSum...')
    s = countTable.sum(axis='columns')
    #plt.ioff()
    #ax = s.plot.hist(bins=20)
    #plot(plt, False, 'stats.png')
    return countTable.div(s, axis='index') * mean(s)

def NormalizeBySpikeIn(countTable):
    print('NormalizeBySpikeIn...')
    s = countTable['Spike1,']
    return countTable.div(s, axis='index')   # TODO: assuming not log values!!!!
    #return countTable.subtract(s, axis='index')   # TODO: assuming log values!!!!


#countTable = pd.DataFrame({'a': [11, 12, 13], 'b': [4, 5, 6]})
#display(countTable)
#display(normalize(countTable))
