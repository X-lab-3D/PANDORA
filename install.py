import os

os.system('mkdir ./PANDORA_files')
os.system('mkdir ./PANDORA_files/data')
os.system('mkdir ./PANDORA_files/data/csv_pkl_files')
os.system('mkdir ./PANDORA_files/outputs')
os.system('mkdir ./PANDORA_files/outputs/logs')
os.system('mkdir ./PANDORA_files/data/PDBs')
os.system('mkdir ./PANDORA_files/data/PDBs/pMHCI')
os.system('mkdir ./PANDORA_files/data/PDBs/IMGT_retrieved/')
os.system('mkdir ./PANDORA_files/data/PDBs/unused_templates')
os.system('mkdir ./PANDORA_files/data/PDBs/unused_templates/parsing_errors')
os.system('mkdir ./PANDORA_files/data/PDBs/unused_templates/non_canonical_res')
os.system('mkdir ./PANDORA_files/data/dist_files')

if "contact-chainID_allAtoms" not in os.listdir('./PANDORA/tools'):
    os.popen('g++ ./PANDORA/tools/contact-chainID_allAtoms.cpp -o ./PANDORA/tools/contact-chainID_allAtoms').read()
else:
    print('contact-chainID_allAtoms script already compiled in PANDORA/tools/')
