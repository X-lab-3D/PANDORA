# 23-Aug-2018 14:36

import modeller as M
import modeller.automodel as MA
from MyLoop import MyLoop
import sys
import os
print()
M.log.verbose()                                # request verbose output
env = M.environ(restyp_lib_file=os.path.join(os.path.dirname(os.getcwd()), "PTM_implementation/combined_restyp_file.lib"))                              # create a new MODELLER environment to build this model in
# directories for input atom files
env.io.atom_files_directory = ['./']

# Read in HETATM records from template PDBs
env.io.hetatm = True

a = MyLoop(env, alnfile= '%s',
              knowns= '%s', sequence = '%s',           #Be sure those two arguments are always in the same line!
              loop_assess_methods = MA.assess.DOPE)
a.toplib = os.path.join(os.path.dirname(os.getcwd()), "PTM_implementation/combined_top_heav_lib_file.lib")
a.make(exit_stage=2)
