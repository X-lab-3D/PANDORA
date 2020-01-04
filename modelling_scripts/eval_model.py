#!/usr/bin/env python
# Li Xue
#  6-Mar-2019 14:49

import sys
import re
from modeller import *
from modeller.scripts import complete_pdb

log.verbose()                                # request verbose output
env = environ()                              # create a new MODELLER environment to build this model in
env.libs.topology.read(file='$(LIB)/top_heav.lib') # read topology
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

# read model file
mdl = complete_pdb(env, 'query_1ogt_AC.BL00070001.pdb')

# Assess with DOPE:
###s = selection(mdl)   # all atom selection
###s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file='query.profile',
###                      normalize_profile=True, smoothing_window=15)

# Select all atoms in the first chain
atmsel = selection(mdl.chains[1])

score = atmsel.assess_dope()
