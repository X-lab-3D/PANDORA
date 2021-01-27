# -*- coding: utf-8 -*-

import os
import sys
from pdb2sql import StructureSimilarity

decoy = sys.argv[1]
ref = sys.argv[2]
#decoy = 'query_1.BL00010001.pdb'
#ref = '1T0M_MP.pdb'
sim = StructureSimilarity(decoy,ref)
#lrmsd_fast = sim.compute_lrmsd_fast(method='svd',
#    lzone='1T0M_MP.lzone')
lrmsd = sim.compute_lrmsd_pdb2sql(method='svd', name = 'CA')

print(lrmsd)
