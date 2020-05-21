import sys
import os

infile = sys.argv[1]
outfile_name = infile + '.temp'
outfile = open(outfile_name, 'w')

for line in open(infile, 'r'):
    s = line.split(' ')
    l = [x for x in s if x != '']
    if len(l) >= 4:
        if l[3] == 'SEP' and l[2] not in ['P', 'O1P', 'O2P', 'O3P', 'HA', 'HB2', 'HB3']:
            outfile.write(line.replace('HETATM', 'ATOM  ').replace('SEP', 'SER'))
        elif l[3] == 'SEP' and l[2] in ['P', 'O1P', 'O2P', 'O3P', 'HA', 'HB2', 'HB3']:
            pass
        else:
            outfile.write(line)
    else:
        outfile.write(line)
outfile.close()

os.system('mv %s %s' %(outfile_name, infile))
