#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 12:22:54 2020

@author: dario
"""

from matplotlib import pyplot as plt
import csv
import sys
import os
from numpy import corrcoef

#indir = sys.argv[1]
os.chdir('../outputs/benchmark')
indir = './'

molpdf= []
DOPE = []
RMSD = []

out = open('RMSD_outliers.txt', 'w')
for fol in os.listdir(indir):
    if os.path.isdir(indir+fol):
        try:
            with open(indir+fol+'/final_scores.tsv') as infile:
                r = csv.reader(infile, delimiter='\t')
                for i, row in enumerate(r):
                    if i == 0:
                        pass
                    else:
                        if float(row[3]) < 4.0:
                            molpdf.append(float(row[1]))
                            DOPE.append(float(row[2]))
                            RMSD.append(float(row[3]))
                        else:
                            out.write(fol + ' ' + row[3] + '\n')
        except:
            pass

out.close()

plt.plot(molpdf, RMSD, 'go')
plt.show()
print('molpdf correlation coefficent:')
print(corrcoef(molpdf, RMSD))

###
# instead corrcoef: 

plt.plot(DOPE, RMSD, 'bo')
plt.show()
print('DOPE correlation coefficent:')
print(corrcoef(DOPE, RMSD))