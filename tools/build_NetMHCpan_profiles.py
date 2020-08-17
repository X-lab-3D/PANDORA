# -*- coding: utf-8 -*-

import os
import re

#%%
'''
os.chdir('../data/allele_profiles')

os.system('wget http://www.cbs.dtu.dk/suppl/immunology/NAR_NetMHCpan_NetMHCIIpan/NetMHCpan_train.tar.gz')
os.system('tar -xzf NetMHCpan_train.tar.gz')
os.system('rm NetMHCpan_train.tar.gz')

os.mkdir('./MSAs')
os.chdir('./NetMHCpan_train/')
'''
#%%

with open('allelelist') as allfile:
    alleles = {}
    for line in allfile:
        row = re.split(' |\t|,', line)
        for allele in row:
            alleles[allele.replace('\n', '')] = []

#print(alleles)

for file in os.listdir('./'):
    if file != 'allelelist' and file != 'MHC_pseudo.dat':
        if file.split('_')[1] == 'ba':
            with open(file) as infile:
                for line in infile:
                    row = line.split(' ')
                    alleles[row[2].replace('\n', '')].append(row[0])
        elif file.split('_')[1] == 'el':
            with open(file) as infile:
                for line in infile:
                    row = line.split(' ')
                    if row[1] == '1':
                        alleles[row[2].replace('\n', '')].append(row[0])
                        
os.chdir('../MSAs/')

for key in alleles:
    with open(key + '.fasta', 'w') as outfile:
        for i, peptide in enumerate(alleles[key]):
            outfile.write('> Seq_' + str(i) + '\n')
            outfile.write(peptide + '\n')
            
                    