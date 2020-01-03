#!/usr/bin/python

import urllib.request
import os
import subprocess
import csv
import gzip
import shutil
import re
import string
import pickle
from Bio.PDB import PDBParser
from Bio.Data.SCOPData import protein_letters_3to1 as to_one_letter_code

### OPEN CSV FILE WITH PDB CODES ###

def imgt_retrieve_clean(ids_filename):

    IDs = []
    with open(ids_filename, 'r') as idsfile:
        spamreader = csv.reader(idsfile, delimiter='\t')
        for i, row in enumerate(spamreader):
            if i == 0:
                pass
            else:
                ID = row[3]
                allele = row[4]
                IDs.append((ID, allele))

    cwd = os.getcwd()
    #print('CWD:', cwd)
    bad_IDs = []
    IDs_dict = {}
    for entry in IDs:
        ID = entry[0]
        ID = ID.upper()
        ID = ID.rstrip()
        url = 'http://www.imgt.org/3Dstructure-DB/IMGT-FILE/IMGT-%s.pdb.gz' %ID
        filepath = '%s/data/PDBs/%s.pdb' %(cwd, ID)
        print('Fetching ', url)

        try:
            urllib.request.urlretrieve( url, filepath + '.gz')
        except: #urllib.error.HTTPError:
            bad_IDs.append(ID)
            print('##################################')
            print('URL not found. Added to Bad_ID')
            print('##################################')
            continue
        with gzip.open('data/PDBs/%s.pdb.gz' %ID, 'rb') as f_in:
            with open('data/PDBs/%s.pdb' %ID, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        subprocess.check_call(['rm', 'data/PDBs/%s.pdb.gz' %ID])
        selected_chains_filepath = '%s/data/PDBs/%s_AC.pdb'%(cwd, ID)
        a_renamed_filepath = '%s/data/PDBs/%s_AA.pdb'%(cwd, ID)
        new_filepath = '%s/data/PDBs/%s_AP.pdb'%(cwd, ID)


        print('Opening %s' %ID)
        pdb = open('data/PDBs/%s.pdb' %ID, 'r')
        aflag = False
        pflag = False
        aID = False
        pID = False
        exflag = False
        count_dict = {}
        alphabet = list(string.ascii_uppercase)
        for line in pdb:
            if not aID:
                if 'COMPND' in line and 'MOLECULE' in line:
                    aflag = True
                    continue
            if not pID:
                if 'COMPND' in line and 'PEPTIDE' in line and not 'PEPTIDE BINDING':
                    pflag = True
                    continue
            if aflag:
                if 'COMPND' in line and 'CHAIN:' in line:
                    newline = line.split(':')[1]
                    aID = re.split("[^a-zA-Z]*", newline)[2]
                    print('aID:', aID)
                    aflag = False
            if pflag:
                if 'COMPND' in line and 'CHAIN:' in line:
                    newline = line.split(':')[1]
                    pID = re.split("[^a-zA-Z]*", newline)[2]
                    print('pID:', pID)
                    pFlag = False
            if line.startswith('ATOM') and 'CA' in line:
                cID = (re.split("[^a-zA-Z]*", line))[13]
                if cID in alphabet:
                    if cID in count_dict:
                        count_dict[cID] += 1
                    else:
                        count_dict[cID] = 1
                else:
                    cID = (re.split("[^a-zA-Z]*", line))[14]
                    if cID in count_dict:
                        count_dict[cID] += 1
                    else:
                        count_dict[cID] = 1

                exflag = True
        if exflag:
            pID = min(count_dict, key=count_dict.get)
            print(count_dict)
            print('pID:', pID)
        count_dict['allele'] = entry[1]
        if ((len(count_dict)-1) %3) == 0:
            IDs_dict[ID] = count_dict
            os.system('pdb_selchain -%s,%s %s > %s' %(aID, pID, filepath, selected_chains_filepath))
            os.system('pdb_rplchain -%s:A %s > %s' %(aID, selected_chains_filepath, a_renamed_filepath))
            os.system('pdb_rplchain -%s:P %s > %s' %(pID, a_renamed_filepath, new_filepath))
        else:
            bad_IDs.append(ID)
        try:
            subprocess.check_call(['rm', 'data/PDBs/%s.pdb' %ID])
        except:
            pass
        try:
            subprocess.check_call(['rm', 'data/PDBs/%s_AC.pdb' %ID])
        except:
            pass
        try:
            subprocess.check_call(['rm', 'data/PDBs/%s_AA.pdb' %ID])
        except:
            pass
        aID = None
        pID = None

    IDd = open("data/IDs_ChainsCounts_dict.pkl", "wb")
    pickle.dump(IDs_dict, IDd)
    IDd.close()
    print('BAD IDs:')
    print(bad_IDs)
    return IDs_dict

def get_pdb_seq(IDs):
    sequences = []
    for ID in IDs:
        seqs = []
        pdbf = 'data/PDBs/%s_AP.pdb' %ID
        P = PDBParser(QUIET=1)
        structure = P.get_structure('s', pdbf)
        for chain in structure.get_chains():
            sequence = ''
            for res in chain:
              try:
                aa = to_one_letter_code[res.resname]
              except KeyError:
                aa= '.'
              sequence += aa
            #print("Chain %s: %s" %(chain.id, sequence))
            seqs.append(sequence)
        sequences.append(seqs)
    return sequences
            
