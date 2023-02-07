import os
from pathlib import Path
from os.path import exists
import json

user_folder_path = Path(__file__).parents[0]

if exists('config.json'):
    with open('config.json') as f:
        data = json.load(f)
        data_folder = data['data_folder_name']
else:
    data_folder = 'default'

dirs = [
        f'{user_folder_path}/Databases', 
        f'{user_folder_path}/Databases/{data_folder}',
        f'{user_folder_path}/Databases/{data_folder}/mhcseqs', 
        f'{user_folder_path}/Databases/{data_folder}/BLAST_databases',
        f'{user_folder_path}/Databases/{data_folder}/PDBs',
        f'{user_folder_path}/Databases/{data_folder}/PDBs/pMHCI', 
        f'{user_folder_path}/Databases/{data_folder}/PDBs/pMHCII',
        f'{user_folder_path}/Databases/{data_folder}/PDBs/Bad', 
        f'{user_folder_path}/Databases/{data_folder}/PDBs/Bad/pMHCI',
        f'{user_folder_path}/Databases/{data_folder}/PDBs/Bad/pMHCII', 
        f'{user_folder_path}/Databases/{data_folder}/PDBs/IMGT_retrieved',
        f'{user_folder_path}/Databases/{data_folder}/outputs',
        f'{user_folder_path}/test/',
        f'{user_folder_path}/test/test_data',
        f'{user_folder_path}/test/test_data/PDBs/Bad',
        f'{user_folder_path}/test/test_data/PDBs/Bad/pMHCI',
        f'{user_folder_path}/test/test_data/PDBs/Bad/pMHCII', 
        ]

for D in dirs:
    try:
        os.mkdir(D)
    except OSError as e:
        print(f'Could not make directory: {D} \n Reason: {e}')

try:
    print('Downloading pre-built database from zenodo...')
    os.popen(f'wget https://sandbox.zenodo.org/record/1129456/files/default.tar.gz?download=1 -O {user_folder_path}/Databases/default.tar.gz').read()
    print('Copying the database')
    os.popen(f'tar -xzvf {user_folder_path}/Databases/default.tar.gz -C {user_folder_path}/Databases/{data_folder}').read()
    os.popen(f'rm {user_folder_path}/Databases/default.tar.gz').read()
except Exception as e:
    print(f'WARNING: received error while installing database: {e}')
    print('To be able to use PANDORA you will have to generate a new database. Please follow the instructions in the README.')