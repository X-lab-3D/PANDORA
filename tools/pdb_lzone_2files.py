#!/usr/bin/env python
# Farzaneh Meimandi Parizi
#  4-Mar-2020 11:57

import sys
import re
from Bio.PDB import PDBParser

ref_pdbf = sys.argv[1]
mob_pdbf = sys.argv[2]
output = open('ref.lzone', 'w')

print(sys.argv[1], sys.argv[2])

P = PDBParser(QUIET=1)
ref_structure = P.get_structure('r', ref_pdbf)
mob_structure = P.get_structure('m', mob_pdbf)
for ref_chain, mob_chain in zip(ref_structure.get_chains(), mob_structure.get_chains()):
    if ref_chain.id == 'M' or ref_chain.id == 'N':
        for ref_residue, mob_residue in zip(ref_chain, mob_chain):
            #print(chain.id, residue.id[1], chain.id, residue.id[1])
            output.write('zone %s%s-%s%s:%s%s-%s%s\n' %(ref_chain.id, str(ref_residue.id[1]),
                ref_chain.id, str(ref_residue.id[1]), mob_chain.id, str(mob_residue.id[1]),
                mob_chain.id, str(mob_residue.id[1])))
    elif ref_chain.id == 'P':
        output.write('fit\n')
        for ref_residue, mob_residue in zip(ref_chain, mob_chain):
            output.write('rzone %s%s-%s%s:%s%s-%s%s\n' %(ref_chain.id, str(ref_residue.id[1]),
                ref_chain.id, str(ref_residue.id[1]), mob_chain.id, str(mob_residue.id[1]),
                mob_chain.id, str(mob_residue.id[1])))
    else:
        raise Exception('Unrecognized chain ID, different from M, N or P. Please check your file')

output.close()
