###  ###
import os
version='1.0'
PANDORA_path = os.path.dirname(os.path.abspath(__file__))
'''
Path to PANDORA installation
:meta hide-value:
'''

PANDORA_data = os.path.join(os.path.dirname(PANDORA_path), 'PANDORA_files', 'data')
'''
Default path to PANDORA data (template database, sequence database, etc.)
:meta hide-value:
'''
