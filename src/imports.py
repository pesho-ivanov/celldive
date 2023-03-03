#from __future__ import print_function
#from IPython.display import display, HTML

#import os
import sys
from pathlib import Path
### MATH ###
import numpy as np
import pandas as pd
#import pylab
import scipy
import scipy.cluster.hierarchy as sch
from scipy.stats import multivariate_normal
import math


# rpy2
#import rpy2.robjects as ro
#import rpy2.robjects.numpy2ri
#ro.numpy2ri.activate()
#R = ro.r


import matplotlib as mpl
#mpl.rcParams['savefig.dpi'] = 300
#mpl.use('Agg')  # or PS, PDF, Agg(for PNG)
mpl.use('SVG')  # or PS, PDF, Agg(for PNG)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#from collections import OrderedDict

#import mpld3; mpld3.enable_notebook()
#import seaborn as sns; #sns.set()

#from pprint import pprint

### PESHO ###
from utils import *
#from file_utils import *
from structs.batch import Batch
from running import Run
from algo.gene_intersect import GeneIntersect
import algos
import test
from samples import *
#from structs.clone_struct import Coloring

if (not "." in sys.path):
    sys.path.append(".")
    
#print sys.path
#print dir(vizualization)

p = print
from IPython.display import display

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
