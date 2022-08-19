###  ###
import os
version='1.0'
PANDORA_path = os.path.dirname(os.path.abspath(__file__))
'''
Path to PANDORA installation
:meta hide-value:
'''

PANDORA_data = os.path.join(os.path.dirname(PANDORA_path), 'Databases', 'default')
# PANDORA_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

alpha_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G',
                'HLA-DQA', 'HLA-DRA', 'HLA-DPA',
                'H2-EA', 'MH2-AA']
beta_genes = ['HLA-DQB', 'HLA-DRB', 'HLA-DPB',
               'H2-EB', 'MH2-AB', 'H2-AB']

MHCI_G_domain=[(1,187)] #182 + 5 for tolerance
MHCII_G_domain=[(1,86),(1,95)] #81 + 5 and 90 + 5 for tolerance