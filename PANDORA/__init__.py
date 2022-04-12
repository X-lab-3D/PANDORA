###  ###
import os
version='0.0.9'
PANDORA_path = os.path.dirname(os.path.abspath(__file__))
PANDORA_data = os.path.join(os.path.dirname(PANDORA_path), 'PANDORA_files', 'data')
# PANDORA_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

alpha_genes = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-E', 'HLA-F', 'HLA-G',
                'HLA-DQA', 'HLA-DRA', 'HLA-DPA', 
                'H2-EA', 'MH2-AA']
beta_genes = ['HLA-DQB', 'HLA-DRB', 'HLA-DPB', 
               'H2-EB', 'MH2-AB', 'H2-AB']
