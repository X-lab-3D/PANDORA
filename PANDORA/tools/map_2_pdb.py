# this script maps and renumbers model.pdb according to ref.pdb
# Here is an example.
import sys
import os

sys.path.append('/home/dariom/PANDORA_package/')
import PANDORA

ref = sys.argv[1]
model= sys.argv[2]

os.system("egrep '^(ATOM|HETATM|END)' %s > clean_%s" %(model, model))

for chnID in ['M','P']:
    os.system(PANDORA.PANDORA_path +'/tools/pdb-pdbalign '+ref+
              ' '+chnID+' clean_'+model+' '+chnID+' > common.pdb')

    # delete the warning line in common.pdb
    os.system("sed -i '/Warning/d' common.pdb")
    
    # delete the lines with residue name of 'X'
    os.system("awk '%s' common.pdb  >  common_%s.pdb" %('$1 == "ATOM" && substr($0,22,1) != "X"', chnID))

os.system('cat common_M.pdb common_P.pdb > common.pdb')
os.system("sed '$ a END' common.pdb")
