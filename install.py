import os

dirs = [
        './PANDORA_files', './PANDORA_files/default', './PANDORA_files/default/csv_pkl_files',
        './PANDORA_files/default/csv_pkl_files/mhcseqs', './PANDORA_files/default/PDBs',
        './PANDORA_files/default/PDBs/pMHCI', './PANDORA_files/default/PDBs/pMHCII',
        './PANDORA_files/default/PDBs/Bad', './PANDORA_files/default/PDBs/Bad/pMHCI',
        './PANDORA_files/default/PDBs/Bad/pMHCII', './PANDORA_files/default/PDBs/IMGT_retrieved',
         './PANDORA_files/default/outputs',
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

