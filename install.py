import os

dirs = [
        './Databases', './Databases/default', './Databases/default/csv_pkl_files',
        './Databases/default/csv_pkl_files/mhcseqs', './Databases/default/PDBs',
        './Databases/default/PDBs/pMHCI', './Databases/default/PDBs/pMHCII',
        './Databases/default/PDBs/Bad', './Databases/default/PDBs/Bad/pMHCI',
        './Databases/default/PDBs/Bad/pMHCII', './Databases/default/PDBs/IMGT_retrieved',
         './Databases/default/outputs',
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

