# -*- coding: utf-8 -*-

import os

filename_start = 'BL00'
filename_end = '0001.pdb'

with open('file.list', 'w') as out:
    for f in os.listdir('.'):
        if filename_start in f and filename_end in f:
            out.write(f + '\n')
            