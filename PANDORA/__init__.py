###  ###
import os
from os.path import exists
from pathlib import Path
import json
# populate package namespace
from PANDORA.Contacts import Contacts
from PANDORA.Database import Database_functions
from PANDORA.Database import Database

version='1.0'
PANDORA_path = os.path.dirname(os.path.abspath(__file__))
user_folder_path = Path(__file__).parents[1]

'''
Path to PANDORA installation
:meta hide-value:
'''

if exists(user_folder_path / 'config.json'):
    with open(user_folder_path / 'config.json') as f:
        data = json.load(f)
        data_folder = data['data_folder_name']
else:
    data_folder = 'default'

PANDORA_data = os.path.join(os.path.dirname(PANDORA_path), 'Databases', data_folder)

alpha_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G',
                'HLA-DQA', 'HLA-DRA', 'HLA-DPA',
                'H2-EA', 'MH2-AA']
beta_genes = ['HLA-DQB', 'HLA-DRB', 'HLA-DPB',
               'H2-EB', 'MH2-AB', 'H2-AB']

MHCI_G_domain=[(1,187)] #182 + 5 for tolerance
MHCII_G_domain=[(1,86),(1,95)] #81 + 5 and 90 + 5 for tolerance