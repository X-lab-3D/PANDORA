# -*- coding: utf-8 -*-

#python

#python retrieve_IMGT_pmh1.py > output.html

import sys

import urllib.parse
import urllib.request


params = { 'ReceptorType' : 'peptide/MH1',
        'type-entry': 'PDB'}


url = "http://www.imgt.org/3Dstructure-DB/cgi/3Dquery.cgi"

data = urllib.parse.urlencode(params)
data = data.encode('ascii') # data should be bytes
req = urllib.request.Request(url, data)

temp_outfile = open('../outputs/test/test_download/test.html', 'w')
with urllib.request.urlopen(req) as response:
    text = response.read().decode('utf-8')
    temp_outfile.write(text)
temp_outfile.close()

temp_outfile = open('../outputs/test/test_download/test.html', 'r')
for line in temp_outfile:
    if 'pdbcode' in line:
        print(line)
temp_outfile.close()
