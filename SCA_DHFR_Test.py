#!/Users/RobbyMiller/opt/anaconda3/bin/python
from __future__ import division

import os
import time
import matplotlib.pyplot as plt
import numpy as np
import copy
import colorsys
import matplotlib.image as mpimg
from IPython.display import display
from IPython.display import Image
import scipy.cluster.hierarchy as sch
from scipy.stats import scoreatpercentile 
import scaTools as sca
import mpld3
import pickle
#import cPickle as pickle
from optparse import OptionParser

if not os.path.exists('Outputs/'): os.makedirs('Outputs/') 


##############################################################
##############################################################
##############################################################


Dseq = list(); Dsca = list(); Dsect = list()
db = pickle.load(open('Outputs/PF00186_full.db','rb'))
Dseq.append(db['sequence'])
Dsca.append(db['sca'])
Dsect.append(db['sector'])
db2 = pickle.load(open('Outputs/DHFR_PEPM3.db', 'rb'))
Dseq.append(db2['sequence'])
Dsca.append(db2['sca'])
Dsect.append(db2['sector'])
N_alg = 2
AlgName = ['PFAM', 'Manual']



##############################################################
##############################################################
##############################################################


ix = 1
plt.rcParams['figure.figsize'] = 9, 15
for k in range(N_alg):
    # List all elements above the diagonal (i<j):
    listS = [Dsca[k]['simMat'][i,j] for i in range(Dsca[k]['simMat'].shape[0]) \
             for j in range(i+1, Dsca[k]['simMat'].shape[1])]
    #Cluster the sequence similarity matrix
    Z = sch.linkage(Dsca[k]['simMat'],method = 'complete', metric = 'cityblock')
    R = sch.dendrogram(Z,no_plot = True)
    ind = map(int, R['ivl'])
    #Plotting
    plt.rcParams['figure.figsize'] = 14, 4 
    plt.subplot(2,2,ix)
    ix += 1
    listS=np.array(listS)
    print([[i] for i in listS[0:10]],Dseq[k]['Npos'])
    plt.hist(x=listS, label=Dseq[k]['Npos']/2)
    plt.xlabel('Pairwise sequence identities', fontsize=14)
    plt.ylabel('Number', fontsize=14)
    plt.subplot(2,2,ix)
    ix += 1
    # print(np.shape(Dsca[k]['simMat'][np.ix_(ind,ind)]))
    plt.show(Dsca[k]['simMat'][np.ix_(ind,ind)], vmin=0, vmax=1)
    plt.colorbar()




