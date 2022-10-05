import os
from os.path import exists
from pathlib import Path
import json

user_folder_path = Path(__file__).parents[1]
print(user_folder_path)

data_folder = 'test'

with open(user_folder_path / 'config.json', 'w') as f:
   f.write('{"data_folder_name" : "%s"}' %data_folder)

# if exists(user_folder_path / 'config.json'):
#     with open(user_folder_path / 'config.json') as f:
#         data = json.load(f)
#         data_folder = data['data_folder_name']
# else:
# data_folder = 'test'

os.system('mkdir ./Databases')
os.system(f'mkdir ./Databases/{data_folder}')
os.system(f'mkdir ./Databases/{data_folder}/PDBs')
#os.system(f'mkdir ./Databases/{data_folder}/PDBs/pMHCI')
#os.system(f'mkdir ./Databases/{data_folder}/PDBs/pMHCII')
os.system(f'mkdir ./Databases/{data_folder}/PDBs/Bad')
os.system(f'mkdir ./Databases/{data_folder}/PDBs/Bad/pMHCI')
os.system(f'mkdir ./Databases/{data_folder}/PDBs/Bad/pMHCII')
os.system(f'mkdir ./Databases/{data_folder}/PDBs/IMGT_retrieved/')

#os.system('mkdir ./test/test_data/PDBs/Bad')
#os.system('mkdir ./test/test_data/PDBs/Bad/pMHCI')
#os.system('mkdir ./test/test_data/PDBs/Bad/pMHCII')
os.system(f'cp -r ./test/test_data/PDBs/* ./Databases/{data_folder}/PDBs/')
os.system(f'cp ./test/test_data/PANDORA_database.pkl ./Databases/{data_folder}/PANDORA_database.pkl')

     


# netMHCII Pan install

# wd = os.path.dirname(os.path.abspath(__file__))
# # Changing working directory
# os.chdir(wd + '/netMHCIIpan-4.0')
# # Downloading data folder for netMHCIIpan
# os.system('wget https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.0/data.tar.gz')
# # Uncompressing
# os.system('gunzip data.tar.gz')
# os.system('tar -xvf data.tar')
# # remove tar file
# os.system('rm data.tar')
# # Make sure the tmp directory is created
# if not os.path.exists('tmp'):
#     os.mkdir('tmp')
#
# # Edit paths in the netMHCIIpan script
# netMHCIIpan_script = []
# with open('netMHCIIpan') as f:
#     for line in f:
#         netMHCIIpan_script.append(line)
# with open('netMHCIIpan', 'w') as f:
#     for line in netMHCIIpan_script:
#         if 'setenv\tNMHOME' in line:
#             f.write('setenv\tNMHOME\t' + wd + '/netMHCIIpan-4.0\n')
#             f.write('setenv\tTMPDIR\t' + wd + '/netMHCIIpan-4.0/tmp\n')
#         else:
#             f.write(line)
#
# # netMHC Pan install
#
# os.chdir(wd + '/netMHCpan-4.1')
# # Downloading data folder for netMHCIIpan
# os.system('wget https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/data.tar.gz ')
# # Uncompressing
# os.system('gunzip data.tar.gz')
# os.system('tar -xvf data.tar')
# # remove tar file
# os.system('rm data.tar')
# # Make sure the tmp directory is created
# if not os.path.exists('tmp'):
#     os.mkdir('tmp')
#
# # Edit paths in the netMHCIIpan script
# netMHCpan_script = []
# with open('netMHCpan') as f:
#     for line in f:
#         netMHCpan_script.append(line)
# with open('netMHCpan', 'w') as f:
#     for line in netMHCpan_script:
#         if 'setenv\tNMHOME' in line:
#             f.write('setenv\tNMHOME\t' + wd + '/netMHCpan-4.1\n')
#             f.write('setenv\tTMPDIR\t' + wd + '/netMHCpan-4.1/tmp\n')
#         else:
#             f.write(line)
#
# # Changing working directory back
# os.chdir(wd)
