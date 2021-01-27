# -*- coding: utf-8 -*-

import sys
import os

indir = sys.argv[1]
filename_start = 'BL00'
filename_end = '0001.pdb'

for directory in os.listdir(indir):
    if os.path.isdir(indir + '/' + directory):
        with open(indir + '/' + directory + '/file.list', 'w') as out:
            for f in os.listdir(indir + '/' + directory):
                if filename_start in f and filename_end in f:
                    out.write(f + '\n')
            