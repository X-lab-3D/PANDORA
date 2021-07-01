import os

try:
    os.mkdir('./PANDORA_files')
    os.mkdir('./PANDORA_files/data')
    os.mkdir('./PANDORA_files/data/csv_pkl_files')
    os.mkdir('./PANDORA_files/data/PDBs')
    os.mkdir('./PANDORA_files/data/PDBs/pMHCI')
    os.mkdir('./PANDORA_files/data/PDBs/pMHCII')
    os.mkdir('./PANDORA_files/data/PDBs/Bad')
    os.mkdir('./PANDORA_files/data/PDBs/Bad/pMHCI')
    os.mkdir('./PANDORA_files/data/PDBs/Bad/pMHCII')
    os.mkdir('./PANDORA_files/data/PDBs/IMGT_retrieved/')
    os.mkdir('./test/test_data/PDBs/Bad')
    os.mkdir('./test/test_data/PDBs/Bad/pMHCI')
    os.mkdir('./test/test_data/PDBs/Bad/pMHCII')
    os.mkdir('./test/test_data/csv_pkl_files')
except OSError:
    pass

'''
#Install dependenciess
os.popen("alias KEY_MODELLER='XXXX'").read()
os.popen("conda install -y -c salilab modeller").read()
os.popen("conda install -y -c conda-forge biopython").read()
os.popen("conda install -y dill").read()
os.popen("conda install -y -c bioconda muscle").read()
os.popen("pip install pdb-tools").read()
os.popen("pip install pdb2sql").read()

os.popen("pip install -e ./").read()
'''
