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
with urllib.request.urlopen(req) as response:
   print(response.read().decode('utf-8'))
