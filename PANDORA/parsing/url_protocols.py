import os
import sys
import urllib.request
import urllib.parse
from pyparsing import nestedExpr
import pickle
from copy import deepcopy

def download_unzip_imgt_structures(del_inn_files = True, del_kabat_files = True):
    '''
    Downloads the complete structural dataset
    from IMGT database: http://www.imgt.org/download/3Dstructure-DB/IMGT3DFlatFiles.tgz
    The only two arguments delete every non-PDB structure file

    Args:
        del_inn_files(bool) : if True (default) deletes all inn files

        del_kabat_files(bool) : if True (default) deletes all kabat files

    '''

    # Changing working directory
    os.chdir('./PANDORA_files/data/PDBs/IMGT_retrieved/')
    # Downloading IMGT dataset
    os.system('wget http://www.imgt.org/download/3Dstructure-DB/IMGT3DFlatFiles.tgz')
    # Uncompressing
    os.system('gunzip IMGT3DFlatFiles.tgz')
    os.system('tar -xvf IMGT3DFlatFiles.tar')

    try:
        os.system('rm IMGT3DFlatFiles.tgz')
    except:
        pass
    os.system('rm IMGT3DFlatFiles.tar')
    # Removing non-PDB files
    if del_inn_files:
        os.system('rm IMGT3DFlatFiles/*.inn.gz')
    if del_kabat_files:
        os.system('rm IMGT3DFlatFiles/*.prot.gz')
    os.chdir('../../../../')

def download_ids_imgt(ReceptorType, out_tsv = False):
    '''
    Querys IMGT with the ReceptorType for PDBs.
    Returns the list of IDs provided by IMGT.

    Args:
        ReceptorType(str) : Receptor query for IMGT

        out_tsv(False or str) : if not False, produces a tsv file names as out_tsv

    Example:
    >>> pMHCI_list = download_ids_imgt('peptide/MH1', 'peptide_MHCI.tsv')

    '''

    #out_tsv = 'pMHCI_IDs_alleles_from_IMGT.tsv'
    params = { 'ReceptorType' : ReceptorType,
             'type-entry': 'PDB'}

    url = "http://www.imgt.org/3Dstructure-DB/cgi/3Dquery.cgi"

    data = urllib.parse.urlencode(params)
    data = data.encode('ascii') # data should be bytes
    req = urllib.request.Request(url, data)

    with urllib.request.urlopen(req) as response:
        text = response.read().decode('utf-8')
        text = text.splitlines()

    IDs_list = []
    IDs_list = [x for x in text if 'href' in x and 'pdbcode' in x]
    IDs_list = [x.split('"') for x in IDs_list]
    IDs_list = [x[3][-4:] for x in IDs_list]

    if out_tsv:
        outfile = open('PANDORA_files/data/csv_pkl_files/' + out_tsv, 'w')
        outfile.write(ReceptorType + ' IMGT IDs\n')
        for ID in IDs_list:
            outfile.write(ID + '\n')
        outfile.close()

    return IDs_list
