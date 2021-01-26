# -*- coding: utf-8 -*-
from Bio.PDB import PDBParser
import os

def move_uncommon_pdbf(indir, outdir, delete_files = False):
    res_list = ['ALA', 'ILE', 'LEU', 'VAL', 'PHE', 'GLY', 'ARG', 'LYS', 'HIS', 'ASN',
                'GLN', 'ASP', 'GLU', 'SER', 'THR', 'TYR', 'MET', 'CYS', 'TRP', 'PRO']
    
    uncommon_pdbs = []
    for pdbf in os.listdir(indir):
        flag = False
        P = PDBParser(QUIET=1)
        try:
            structure = P.get_structure('s', indir + pdbf)
        except:
            print('Something wrong in the parsing. Check ' + pdbf)
            continue
        for chain in structure.get_chains():
            if chain.id == 'P':
                for i, res in enumerate(chain):
                    if res.resname not in res_list:
                        print(res.resname, chain.id, pdbf)
                        uncommon_pdbs.append(pdbf[0:4])
                        flag = True
                        break
            if flag:
                break
        if flag:
            if delete_files:
                os.system('rm %s/%s' %(indir, pdbf))
            else:
                os.system('mv %s/%s %s' %(indir, pdbf, outdir))
    
    return list(set(uncommon_pdbs))