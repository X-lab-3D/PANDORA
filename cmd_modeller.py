#!/usr/bin/env python
# 23-Aug-2018 14:36

from modeller import *
from modeller.automodel import *
from MyLoop import MyLoop

log.verbose()                                # request verbose output
env = environ()                              # create a new MODELLER environment to build this model in

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']

# Read in HETATM records from template PDBs
#env.io.hetatm = True
env.io.water = True

a = MyLoop(env, alnfile='1k5n_1ogt.ali',
              knowns='1k5n', sequence='query_1ogt_AC',
              loop_assess_methods=assess.DOPE)


a.starting_model= 1                          # index of the first model
a.ending_model  = 2                         # index of the last model

# Optimization
# CG
a.library_schedule = autosched.slow           # Very thorough VTFM optimization
a.max_var_iterations = 300                    # Select length of optimizations
a.max_molpdf = 1e6                            # do not stop unless obj.func. > 1E6

# Loop Modelling

a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 10          # Last loop model
a.loop.md_level       = refine.slow # Loop model refinement level

#MD
a.md_level = refine.slow                      # model refinement level

# Repeat the whole cycle 2 times
#a.repeat_optimization = 2

a.make()                                     # do the actual homology modeling
