# 23-Aug-2018 14:36

import modeller as M
import modeller.automodel as MA
from MyLoop import MyLoop
import sys

alifile = sys.argv[1]
ID = sys.argv[2]
query = sys.argv[3]
#anchor_1 = sys.argv[4]
#anchor_2 = sys.argv[5]
#alifile, ID, query = 'data/Alignments/4PGB.ali', '4PGB', '>3ROO:A'

M.log.verbose()                                # request verbose output
env = M.environ()                              # create a new MODELLER environment to build this model in

# directories for input atom files
#env.io.atom_files_directory = ['.', '../atom_files']
env.io.atom_files_directory = './'

# Read in HETATM records from template PDBs
#env.io.hetatm = True
env.io.water = True

a = MyLoop(env, alnfile=alifile,
              knowns=ID, sequence = query,
              loop_assess_methods = MA.assess.DOPE)


a.starting_model= 1                          # index of the first model
a.ending_model  = 1                         # index of the last model

# Optimization
# CG
a.library_schedule = MA.autosched.slow           # Very thorough VTFM optimization
a.max_var_iterations = 300                    # Select length of optimizations
a.max_molpdf = 1e6                            # do not stop unless obj.func. > 1E6

# Loop Modelling

a.loop.starting_model = 1           # First loop model
a.loop.ending_model   = 10          # Last loop model
a.loop.md_level       = MA.refine.slow # Loop model refinement level

#MD
a.md_level = MA.refine.slow                      # model refinement level

# Repeat the whole cycle 2 times
#a.repeat_optimization = 2

a.make()                                     # do the actual homology modeling