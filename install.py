import os

dirs = [
        './PANDORA_files', './PANDORA_files/data', './PANDORA_files/data/csv_pkl_files', 
        './PANDORA_files/data/csv_pkl_files/mhcseqs', './PANDORA_files/data/PDBs', 
        './PANDORA_files/data/PDBs/pMHCI', './PANDORA_files/data/PDBs/pMHCII', 
        './PANDORA_files/data/PDBs/Bad', './PANDORA_files/data/PDBs/Bad/pMHCI',
        './PANDORA_files/data/PDBs/Bad/pMHCII', './PANDORA_files/data/PDBs/IMGT_retrieved',
        './test/test_data/PDBs/Bad', './test/test_data/PDBs/Bad/pMHCI', 
        './test/test_data/PDBs/Bad/pMHCII', './test/test_data/csv_pkl_files',
        ]
    
for D in dirs:
    try:
        os.mkdir(D)
    except OSError:
        print('Could not make directory: ' + D)

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
