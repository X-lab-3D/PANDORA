#!/usr/bin/env python
# Farzaneh Meimandi Parizi
#  4-Mar-2020 11:57

import sys
from Bio.PDB import PDBParser

ref_pdbf = sys.argv[1]
mob_pdbf = sys.argv[2]
output = open('ref.lzone', 'w')

print(sys.argv[1], sys.argv[2])

P = PDBParser(QUIET=1)
ref_structure = P.get_structure('r', ref_pdbf)
mob_structure = P.get_structure('m', mob_pdbf)
r_str_dict = {}
m_str_dict = {}
for ref_chain in ref_structure.get_chains():
    if ref_chain.id == 'M' or ref_chain.id == 'N' or ref_chain.id == 'P': 
        r_str_dict[ref_chain.id] = [x.id[1] for x in ref_chain if x[2] == ' ']
    
for mob_chain in mob_structure.get_chains():
    if mob_chain.id == 'M' or mob_chain.id == 'N' or mob_chain.id == 'P': 
        m_str_dict[mob_chain.id] = [x.id[1] for x in mob_chain if x[2] == ' ']
    
for key in r_str_dict.keys():
    if key == 'M' or key == 'N':
        for ref_residue, mob_residue in zip(r_str_dict[key], m_str_dict[key]):
            #print(chain.id, residue.id[1], chain.id, residue.id[1])
            output.write('zone %s%s-%s%s:%s%s-%s%s\n' %(key, str(r_str_dict[key]),
                key, str(r_str_dict[key]), key, str(m_str_dict[key]),
                key, str(m_str_dict[key])))
    elif key == 'P':
        output.write('fit\n')
        for ref_residue, mob_residue in zip(r_str_dict[key], m_str_dict[key]):
            output.write('rzone %s%s-%s%s:%s%s-%s%s\n' %(key, str(r_str_dict[key]),
                key, str(r_str_dict[key]), key, str(m_str_dict[key]),
                key, str(m_str_dict[key])))
    else:
        raise Exception('Unrecognized chain ID, different from M, N or P. Please check your file')

output.close()
