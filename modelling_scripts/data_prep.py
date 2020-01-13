#!/usr/bin/python

import urllib.request
import os
import subprocess
import csv
import gzip
import shutil
import pickle
from Bio.PDB import PDBParser
from Bio.Data.SCOPData import protein_letters_3to1 as to_one_letter_code

### OPEN CSV FILE WITH PDB CODES ###

def imgt_retrieve_clean(ids_filename, id_clmn, allele_clmn, delimiter):

    IDs = []
    with open(ids_filename, 'r') as idsfile:
        spamreader = csv.reader(idsfile, delimiter=delimiter)
        for i, row in enumerate(spamreader):
            if i == 0:
                pass
            else:
                ID = row[id_clmn]
                allele = row[allele_clmn]
                IDs.append((ID, allele))

    cwd = os.getcwd()
    #print('CWD:', cwd)
    bad_IDs = []
    IDs_dict = {}
    print('Started fetching URLs')
    for entry in IDs:
        ID = entry[0]
        ID = ID.upper()
        ID = ID.rstrip()
        url = 'http://www.imgt.org/3Dstructure-DB/IMGT-FILE/IMGT-%s.pdb.gz' %ID
        filepath = '%s/data/PDBs/%s.pdb' %(cwd, ID)
        #print('Fetching ', url)

        try:
            urllib.request.urlretrieve( url, filepath + '.gz')
        except: #urllib.error.HTTPError:
            bad_IDs.append(ID)
            print('URL not found. %s added to Bad_IDs' %ID)
            print('##################################')
            continue
        with gzip.open('data/PDBs/%s.pdb.gz' %ID, 'rb') as f_in:
            with open('data/PDBs/%s.pdb' %ID, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        subprocess.check_call(['rm', 'data/PDBs/%s.pdb.gz' %ID])
        selected_chains_filepath = '%s/data/PDBs/%s_AC.pdb'%(cwd, ID)
        a_renamed_filepath = '%s/data/PDBs/%s_MC.pdb'%(cwd, ID)
        new_filepath = '%s/data/PDBs/%s_MP.pdb'%(cwd, ID)


        print('Opening %s' %ID)
        #pdb = open('data/PDBs/%s.pdb' %ID, 'r')
        pdbf = 'data/PDBs/%s.pdb' %ID
        #aflag = False
        #pflag = False
        aID = False
        pID = False
        #exflag = False
        count_dict = {}
        #alphabet = list(string.ascii_uppercase)
        
        try:
            seqs = get_seqs(pdbf)
        except IndexError:
            bad_IDs.append(ID)
            print('Something went wrong with PDBParser. %s added to Bad_IDs' %ID)
            print('##################################')
            subprocess.check_call(['rm', 'data/PDBs/%s.pdb' %ID])
            continue
        for chain in seqs:
            length = len(seqs[chain])
            if chain != ' ':
                count_dict[chain] = length
            if length > 250 and length < 300 and not aID:
                aID = chain
                #print('aID : ', aID)
                #print('length : ', length)
            elif length > 7 and length < 25 and not pID:
                pID = chain
                #print('pID : ', pID)
                #print('length : ', length)
        try:
            if aID.islower() or pID.islower():
                bad_IDs.append(ID)
                print('Lower case chain ID. %s added to Bad_IDs' %ID)
                print('##################################')
                subprocess.check_call(['rm', 'data/PDBs/%s.pdb' %ID])
                continue
        except AttributeError:
            bad_IDs.append(ID)
            print('False in chain ID, one chain does not respect settled lenght criteria. %s added to Bad_IDs' %ID)
            print('##################################')
            subprocess.check_call(['rm', 'data/PDBs/%s.pdb' %ID])
            continue

        count_dict['allele'] = entry[1]
        if ((len(count_dict)-1) %3) == 0:
            IDs_dict[ID] = count_dict
            os.system('pdb_selchain -%s,%s %s > %s' %(aID, pID, filepath, selected_chains_filepath))
            if pID == 'M':
                os.system('pdb_rplchain -%s:P %s > %s' %(pID, selected_chains_filepath, a_renamed_filepath))
                os.system('pdb_rplchain -%s:M %s > %s' %(aID, a_renamed_filepath, new_filepath))
            else:
                os.system('pdb_rplchain -%s:M %s > %s' %(aID, selected_chains_filepath, a_renamed_filepath))
                os.system('pdb_rplchain -%s:P %s > %s' %(pID, a_renamed_filepath, new_filepath))
        else:
            bad_IDs.append(ID)
            print('Chain number not multiple of three. %s added to Bad_IDs' %ID)
            print('##################################')
            continue
        try:
            subprocess.check_call(['rm', 'data/PDBs/%s.pdb' %ID])
        except:
            pass
        try:
            subprocess.check_call(['rm', 'data/PDBs/%s_AC.pdb' %ID])
        except:
            pass
        try:
            subprocess.check_call(['rm', 'data/PDBs/%s_MC.pdb' %ID])
        except:
            pass

    IDd = open("data/IDs_ChainsCounts_dict.pkl", "wb")
    pickle.dump(IDs_dict, IDd)
    pickle.dump(bad_IDs, IDd)
    IDd.close()
    print('BAD IDs:')
    print(bad_IDs)
    with open('data/IDs_dict.tsv', 'wt') as outfile:
        tsv_writer = csv.writer(outfile, delimiter='\t')
        tsv_writer.writerow(['PDB_ID', 'CHAIN_ID', 'LENGTH', 'ALLELE'])
        for ID in IDs_dict:
            for chain_ID in IDs_dict[ID]:
                if chain_ID != 'allele':
                    tsv_writer.writerow([ID, chain_ID, IDs_dict[ID][chain_ID], IDs_dict[ID]['allele']])
    
    return IDs_dict, bad_IDs

def get_pdb_seq(IDs):
    sequences = []
    empty_seqs = []
    for ID in IDs:
        pdbf = 'data/PDBs/%s_MP.pdb' %ID
        seqs = get_seqs(pdbf)
        if seqs != {}:
            sequences.append(seqs)
        else:
            empty_seqs.append(ID)
    #return sequences1, sequences2, empty_seqs
    return sequences, empty_seqs

def get_seqs(pdbf):
    seqs = {}
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
        seqs[chain.id] = sequence
    return seqs
            