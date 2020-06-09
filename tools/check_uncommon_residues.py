# -*- coding: utf-8 -*-
from Bio.PDB import PDBParser
import os

res_list = ['ALA', 'ILE', 'LEU', 'VAL', 'PHE', 'GLY', 'ARG', 'LYS', 'HIS', 'ASN',
            'GLN', 'ASP', 'GLU', 'SER', 'THR', 'TYR', 'MET', 'CYS', 'TRP', 'PRO']

rr = {}
for pdbf in os.listdir('../data/PDBs/'):
    P = PDBParser(QUIET=1)
    try:
        structure = P.get_structure('s', '../data/PDBs/' + pdbf)
    except:
        #print('Something wrong in the parsing. Check ' + pdbf)
        continue
    for chain in structure.get_chains():
        if chain.id != ' ':
            for i, res in enumerate(chain):
                if res.resname not in res_list:
                    print(res.resname, chain.id, pdbf)
                    if chain.id == 'P':
                        try:
                            rr[res.resname] += 1
                        except:
                            rr[res.resname] = 1
