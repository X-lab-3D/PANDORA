from os.path import exists
import os
import json

if exists('config.json'):
    with open('config.json') as f:
        data = json.load(f)
        data_folder = data['data_folder_name']
else:
    data_folder = 'default'

dirs = [
        './Databases', f'./Databases/{data_folder}',
        f'./Databases/{data_folder}/mhcseqs', f'./Databases/{data_folder}/PDBs',
        f'./Databases/{data_folder}/PDBs/pMHCI', f'./Databases/{data_folder}/PDBs/pMHCII',
        f'./Databases/{data_folder}/PDBs/Bad', f'./Databases/{data_folder}/PDBs/Bad/pMHCI',
        f'./Databases/{data_folder}/PDBs/Bad/pMHCII', f'./Databases/{data_folder}/PDBs/IMGT_retrieved',
        f'./Databases/{data_folder}/outputs',
        './test/test_data/PDBs/Bad','./test/test_data/PDBs/Bad/pMHCI',
        './test/test_data/PDBs/Bad/pMHCII', './test/test_data'
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

