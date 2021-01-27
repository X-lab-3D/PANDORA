import os
import sys

folder = sys.argv[1]
filename_start = 'BL00'
filename_end = '0001.pdb'

for f in os.listdir(folder):
    flag = False
    counter = 0
    if '_query_' in f:
        for name in os.listdir('%s/%s' %(folder, f)):
            if filename_start in name and filename_end in name:
                counter += 1
        if counter == 20:
           flag = True
        elif counter == 0:
           print(f, counter)
           os.system('rm -r %s/%s' %(folder, f))
