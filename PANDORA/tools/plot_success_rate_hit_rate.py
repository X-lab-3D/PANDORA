from matplotlib import pyplot as plt
import csv
import sys
import os
import torch
from torch.nn.functional import softmax
from numpy import array
import math
from statistics import mean


####
# LOGIC: input csv with target (0/1) column, first scoring method column, second scoring method column.
# output: success rate curve plot


#Meant to be runned from inside the output/benchmark<mybenchmark> directory

indir = sys.argv[1]
outdir = sys.argv[2]
#outdir = '../../benchmark/no_prize/1906_unlabeled/'

res_dict = {}
for fol in os.listdir(indir):
    if os.path.isdir(indir + fol):
        case = fol.split('_')[1]
        try:
            infile = open(indir + fol + '/' + 'rmsds_and_final_scores.tsv')
        except:
            print(indir + fol + '/' + 'rmsds_and_final_scores.tsv')
            continue
        res_dict[case] = {}
        for i, line in enumerate(infile):
            row = line.split('\t')
            if i == 0:
                headers = row
            else:
                if len(row) == 6:
                    res_dict[case][row[0]] = [float(x.replace('\n', '')) for x in row[1:]]
                else:
                    pass
        if len(res_dict[case]) == 19:
            decoy_model = []
            for value in range(len(row)-1):
                decoy_model.append(mean([res_dict[case][model][value] for model in res_dict[case]]))
                #[[mean(res_dict[case][model][j]) for j in range(len(res_dict[case][model]))] for model in res_dict[case]]
            res_dict[case]['decoy_model'] = decoy_model
        infile.close()


outfile = open(outdir + 'hitrate_rawdata.tsv', 'w')
spamwriter = csv.writer(outfile, delimiter='\t')
spamwriter.writerow(['caseID', 'modelID', 'RMSD', 'target', 'DR', 'HS'])
for case in res_dict:
    for model in res_dict[case]:
        target = 0
        if res_dict[case][model][4] <= 2:
            target = 1
        spamwriter.writerow([case, model, res_dict[case][model][4], target, res_dict[case][model][0], res_dict[case][model][1]])
outfile.close()

#%%
res_tens_dict = {}
for key in res_dict:
    res_tens_dict[key] = torch.Tensor([res_dict[key][x] for x in res_dict[key]])

#res_tens = torch.Tensor(1, 20, 5)
res_tens = None
for x in res_tens_dict:
    sqtens = torch.unsqueeze(res_tens_dict[x], 0)
    try:
        res_tens = torch.cat((res_tens, sqtens), 0)
    except:
        res_tens = sqtens.clone()

res_tens[:,:,0] = torch.sigmoid((res_tens[:,:,0] - res_tens[:,:,0].mean(1).reshape(-1,1))/(res_tens[:,:,0].std(1).reshape(-1,1)))
res_tens[:,:,1] = torch.sigmoid((res_tens[:,:,1] - res_tens[:,:,1].mean(1).reshape(-1,1))/(res_tens[:,:,1].std(1).reshape(-1,1)))

#plt.hist(torch.clamp(res_tens[:,:,0].reshape(-1),-600,75), bins=100)
plt.hist(res_tens[:,:,1].reshape(-1), bins=100)
#plt.savefig(outdir + 'test.pdf')
#%%
'''
for key in res_dict:
    molpdfs = torch.Tensor([res_dict[key][x][0] for x in res_dict[key]])
    dopes = torch.Tensor([res_dict[key][x][1] for x in res_dict[key]])
    norm_molpdfs = 1 / (1 + torch.exp((molpdfs - torch.mean(molpdfs)) / torch.std(molpdfs)))
    norm_dopes = 1 / (1 + torch.exp((dopes - torch.mean(dopes)) / torch.std(dopes)))
    break

#tensor[feature_position] = torch.sigmoid(tensor[feature_position])
    
plt.plot(norm_molpdfs, norm_dopes, 'ro')
'''
