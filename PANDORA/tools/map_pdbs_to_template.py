### USAGE: python map_pdbs_to_template.py <ref_PDB.pdb> <PDBs_list.list/txt>

import sys

template_PDB = sys.argv[1]
PDBs_list = sys.argv[2]

for line in open(PDBs_list):
    f = line.replace('/n', '')
    if f.endswith('.pdb'):
        os.popen('bash ./map_2_pdb.sh %s %s >  matched_%s' %(template_PDB, f,f)).read() #Check for the path to map_2_pdb.sh script!
