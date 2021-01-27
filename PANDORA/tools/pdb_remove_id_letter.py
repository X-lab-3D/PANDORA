# -*- coding: utf-8 -*-

import sys

infile = sys.argv[1]
for line in open(infile, 'r'):
    print((line[:26] + ' ' + line[26 + 1:]).rstrip())