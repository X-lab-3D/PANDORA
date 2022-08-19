import os

dirs = [
        './Databases', './Databases/data', './Databases/data/csv_pkl_files',
        './Databases/data/csv_pkl_files/mhcseqs', './Databases/data/PDBs',
        './Databases/data/PDBs/pMHCI', './Databases/data/PDBs/pMHCII',
        './Databases/data/PDBs/Bad', './Databases/data/PDBs/Bad/pMHCI',
        './Databases/data/PDBs/Bad/pMHCII', './Databases/data/PDBs/IMGT_retrieved',
         './Databases/data/outputs',
        './test/test_data/PDBs/Bad','./test/test_data/PDBs/Bad/pMHCI',
        './test/test_data/PDBs/Bad/pMHCII', './test/test_data/csv_pkl_files'
        ]

for D in dirs:
    try:
        os.mkdir(D)
    except OSError:
        print('Could not make directory: ' + D)


# Install dependenciess
# os.popen("alias KEY_MODELLER='XXXX'").read()
# os.popen("conda install -y -c salilab modeller").read()
# os.popen("conda install -y -c bioconda muscle").read()

# os.popen("pip install -e ./").read()

