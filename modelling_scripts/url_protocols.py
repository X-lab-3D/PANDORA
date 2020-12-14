import os
import sys
import urllib.request
import urllib.parse
from pyparsing import nestedExpr
import pickle
from copy import deepcopy

def download_unzip_imgt_structures(del_inn_files = True, del_kabat_files = True):
    # Changing working directory
    os.chdir('./data/PDBs/IMGT_retrieved/')
    # Downloading IMGT dataset
    os.system('wget http://www.imgt.org/download/3Dstructure-DB/IMGT3DFlatFiles.tgz')
    # Uncompressing
    os.system('gunzip IMGT3DFlatFiles.tgz')
    os.system('tar -xvf IMGT3DFlatFiles.tar')
    
    os.system('rm IMGT3DFlatFiles.tar')
    # Removing non-PDB files
    if del_inn_files:
        os.system('rm IMGT3DFlatFiles/*.inn.gz')
    if del_kabat_files:
        os.system('rm IMGT3DFlatFiles/*.prot.gz')

def download_ids_imgt(ReceptorType, out_tsv = 'auto_generated_IDs_alleles_from_IMGT.tsv', out_pkl = 'IDs_and_alleles_identity_percs_from_imgt.pkl', print_outfiles = False):
    
    '''
    params = { 'ReceptorType' : 'peptide/MH1',
            'type-entry': 'PDB'}
    '''
    params = { 'ReceptorType' : ReceptorType,
        'type-entry': 'PDB'}
    
    url = "http://www.imgt.org/3Dstructure-DB/cgi/3Dquery.cgi"
    
    data = urllib.parse.urlencode(params)
    data = data.encode('ascii') # data should be bytes
    req = urllib.request.Request(url, data)
    
    with urllib.request.urlopen(req) as response:
        text = response.read().decode('utf-8')
        text = text.splitlines()
        if print_outfiles: 
            temp_outfile = open('./outputs/test/test_download/test.html', 'w')
            temp_outfile.write(text)
            temp_outfile.close()
    
    IDs_list = []
    IDs_list = [x for x in text if 'href' in x and 'pdbcode' in x]
    IDs_list = [x.split('"') for x in IDs_list]
    IDs_list = [x[3][-4:] for x in IDs_list]
    
    return IDs_list