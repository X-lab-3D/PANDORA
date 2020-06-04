# -*- coding: utf-8 -*-

import os
import sys

clmn = sys.argv[1]

#print(clmns)
bests = []
for fol in os.listdir('./'):
    #print(fol)
    if os.path.isdir(fol):
        if os.path.isfile(fol + '/rmsds_and_final_scores.tsv'):
            rmsds = []
            with open(fol + '/rmsds_and_final_scores.tsv')as infile:
                for i, line in enumerate(infile):
                    if i != 0:
                        rmsds.append(float(line.split('\t')[int(clmn)]))
                        #print(float(line.split('\t')[int(col)]))
            best = sorted(rmsds)[:5]
            #print('BEST: ', sum(best)/len(best))
            bests.append(sum(best)/len(best))

print ('Final AVG: ', sum(bests)/len(bests))