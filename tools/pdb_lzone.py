#!/usr/bin/env python
# Farzaneh Meimandi Parizi
#  4-Mar-2020 11:57

import sys
import re
from Bio.PDB import PDBParser

ref_pdbf = sys.argv[1]
output = open('ref.lzone', 'w')

P = PDBParser(QUIET=1)
structure = P.get_structure('r', ref_pdbf)
for chain in structure.get_chains():
    if chain.id == 'M' or chain.id == 'N':
        for residue in chain:
            output.write('zone %s%s-%s%s\n' %(chain.id, str(residue.id[1]), chain.id, str(residue.id[1])))
    elif chain.id == 'P':
        output.write('fit\n')
        for residue in chain:
            output.write('rzone %s%s-%s%s\n' %(chain.id, str(residue.id[1]), chain.id, str(residue.id[1])))
    else:
        raise Exception('Unrecognized chain ID, different from M, N or P. Please check your file')

output.close()
