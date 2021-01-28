#!/usr/bin/python

import urllib.request
import urllib.parse
from pyparsing import nestedExpr
from copy import deepcopy
import os
import subprocess
import csv
import gzip
import shutil
import pickle

from numpy import argmax
from Bio.PDB import PDBParser
from Bio.Data.SCOPData import protein_letters_3to1 as to_one_letter_code


'''
This file contains the following functions:
    - align_peptides(seq1, anch1_seq1, anch2_seq1, seq2, anch1_seq2, anch2_seq2)
    - get_peptides_from_csv(pepts_filename, pept_clmn, allele_clmn, delimiter, skip_first_line = True)
    - get_peptides_w_star_from_csv(pepts_filename, pept_clmn, allele_clmn, delimiter, skip_first_line = True)
    - remove_error_templates(pklfile)
    - change_sep_in_ser(indir)
    - get_pdb_seq(IDs)
    - get_seqs(pdbf)
'''


### OPEN CSV FILE WITH PDB CODES ###
#def pdb_reres_allchains(pdb_file)

def align_peptides(seq1, anch1_seq1, anch2_seq1, seq2, anch1_seq2, anch2_seq2):
    '''
    Align two MHC-I peptides making overlap the anchors.
    This function does NOT use an alignment matrix (e.g. BLOSUM, PAM, etc).
    It computes a simple anchor position alignment and inserts gap in the
    middle part to make the final sequences have the same lenghts.

    Args:
        seq1(str) : sequence of the first peptide.
        anch1_seq1(int) : position of the first anchor of seq1. Position must be given in Python numbering (0-N)
        anch2_seq1(int) : position of the second anchor of seq1. Position must be given in Python numbering (0-N)
        seq2(str) : sequence of the second peptide.
        anch1_seq1(int) : position of the first anchor of seq1. Position must be given in Python numbering (0-N)
        anch2_seq1(int) : position of the second anchor of seq1. Position must be given in Python numbering (0-N)

    Returns:
        ali_seq1(str)
    '''

    seq1_core = anch2_seq1 - anch1_seq1
    seq2_core = anch2_seq2 - anch1_seq2
    tail1 = [x for x in seq1[anch2_seq1:]]
    tail2 = [x for x in seq2[anch1_seq2:]]

    list1 = [x for x in seq1]
    list2 = [x for x in seq2]
    #Adding gaps in cores
    if seq1_core > seq2_core:
        for x in range(seq1_core - seq2_core):
            list2.insert(int(len(seq2)/2), '-')
    elif seq1_core < seq2_core:
        for x in range(seq2_core - seq1_core):
            list1.insert(int(len(seq1)/2), '-')
    ### Adding gaps in heads
    if anch1_seq1 > anch1_seq2:
        for x in range(anch1_seq1 - anch1_seq2):
            list2.insert(0, '-')
    elif anch1_seq1 < anch1_seq2:
        for x in range(anch1_seq2 - anch1_seq1):
            list1.insert(0, '-')
    ### Adding gaps in heads
    if len(tail1) > len(tail2):
        for x in range(len(tail1) - len(tail2)):
            list2.insert(-1, '-')
    elif len(tail1) < len(tail2):
        for x in range(len(tail1) - len(tail2)):
            list1.insert(-1, '-')

    ali_seq1 = ('').join(list1)
    ali_seq2 = ('').join(list2)
    return ali_seq1, ali_seq2

def get_peptides_from_csv(pepts_filename, pept_clmn, allele_clmn, delimiter, skip_first_line = True):

    seqs = []
    with open(pepts_filename, 'r') as peptsfile:
        spamreader = csv.reader(peptsfile, delimiter=delimiter)
        for i, row in enumerate(spamreader):
            if i == 0 and skip_first_line:
                pass
            else:
                seq = row[pept_clmn]
                allele = row[allele_clmn]
                seqs.append((seq, allele))
    return seqs

def get_peptides_w_star_from_csv(pepts_filename, pept_clmn, allele_clmn, delimiter, skip_first_line = True):

    seqs = []
    with open(pepts_filename, 'r') as peptsfile:
        spamreader = csv.reader(peptsfile, delimiter=delimiter)
        for i, row in enumerate(spamreader):
            if i == 0 and skip_first_line:
                pass
            else:
                seq = row[pept_clmn]
                allele = row[allele_clmn]
                if 'HLA' in allele:
                    star_allele = (allele[0:5]+'*'+allele[5:])
                    seqs.append((seq, star_allele))
                else:
                    seqs.append((seq, allele))
    return seqs

def remove_error_templates(pklfile):
    try:
        os.system('mkdir PANDORA_files/data/PDBs/unused_templates')
    except:
        pass
    with open('PANDORA_files/data/csv_pkl_files/%s' %pklfile, 'rb') as inpkl:
        IDD = pickle.load(inpkl)
        bad_IDs = pickle.load(inpkl)
        inpkl.close()
    with open('PANDORA_files/data/csv_pkl_files/error_templates.tsv') as infile: #To be removed later
        r = csv.reader(infile, delimiter='\t')
        for i, row in enumerate(r):
            if i != 0:
                try:
                    del IDD[row[0]]
                    bad_IDs[row[0]] = row[1]
                    os.system('mv PANDORA_files/data/PDBs/pMHCI/%s_MP.pdb PANDORA_files/data/PDBs/unused_templates/' %row[0])
                except:
                    pass
        infile.close()
    with open("PANDORA_files/data/csv_pkl_files/%s" %pklfile, "wb") as outpkl:
        pickle.dump(IDD, outpkl)
        pickle.dump(bad_IDs, outpkl)
        outpkl.close()
        return IDD, bad_IDs

def change_sep_in_ser(indir):
    for pdbfile in os.listdir(indir):
        outfile_name = indir + pdbfile + '.temp'
        outfile = open(outfile_name, 'w')
        for line in open(indir + '/' + pdbfile, 'r'):
            s = line.split(' ')
            l = [x for x in s if x != '']
            if line.startswith('ATOM') or line.startswith('HETATM'):
                if ('SEP' in l[3] or 'SEP' in l[2]) and l[2] not in ['P', 'O1P', 'O2P', 'O3P', 'HA', 'HB2', 'HB3']: #Lines to keep
                    outfile.write(line.replace('HETATM', 'ATOM  ').replace('SEP', 'SER'))
                elif ('SEP' in l[3] or 'SEP' in l[2]) and l[2] in ['P', 'O1P', 'O2P', 'O3P', 'HA', 'HB2', 'HB3']:   #Lines to delete
                    pass
                elif ('F2F' in l[3] or 'F2F' in l[2]) and l[2] not in ['F1', 'F2']:
                    outfile.write(line.replace('HETATM', 'ATOM  ').replace('F2F', 'PHE'))
                elif ('F2F' in l[3] or 'F2F' in l[2]) and l[2] in ['F1', 'F2']:
                    pass
                elif ('CSO' in l[3] or 'CSO' in l[2]) and l[2] not in ['OD']: #Lines to keep
                    outfile.write(line.replace('HETATM', 'ATOM  ').replace('CSO', 'CYS'))
                elif ('CSO' in l[3] or 'CSO' in l[2]) and l[2] in ['OD']:   #Lines to delete
                    pass
                else:
                    outfile.write(line)
            else:
                outfile.write(line)
        outfile.close()

        os.system('mv %s %s' %(outfile_name, indir + pdbfile))

def move_uncommon_pdbf(indir, outdir, delete_files = False):
    res_list = ['ALA', 'ILE', 'LEU', 'VAL', 'PHE', 'GLY', 'ARG', 'LYS', 'HIS', 'ASN',
                'GLN', 'ASP', 'GLU', 'SER', 'THR', 'TYR', 'MET', 'CYS', 'TRP', 'PRO']

    uncommon_pdbs = []
    for pdbf in os.listdir(indir):
        flag = False
        P = PDBParser(QUIET=1)
        try:
            structure = P.get_structure('s', indir + pdbf)
        except:
            print('Something wrong in the parsing. Check ' + pdbf)
            continue
        for chain in structure.get_chains():
            if chain.id == 'P':
                for i, res in enumerate(chain):
                    if res.resname not in res_list:
                        print(res.resname, chain.id, pdbf)
                        uncommon_pdbs.append(pdbf[0:4])
                        flag = True
                        break
            if flag:
                break
        if flag:
            if delete_files:
                os.system('rm %s/%s' %(indir, pdbf))
            else:
                os.system('mv %s/%s %s/' %(indir, pdbf, outdir))

    return list(set(uncommon_pdbs))

def get_pdb_seq(IDs):
    sequences = []
    empty_seqs = []
    for ID in IDs:
        pdbf = 'PANDORA_files/data/PDBs/pMHCI/%s_MP.pdb' %ID
        seqs = get_seqs(pdbf)
        if seqs != {}:
            sequences.append(seqs)
        else:
            empty_seqs.append(ID)
    return sequences, empty_seqs

def get_seqs(pdbf):
    seqs = {}
    P = PDBParser(QUIET=1)
    structure = P.get_structure('s', pdbf)
    for chain in structure.get_chains():
        if chain.id == ' ':
            continue
        sequence = ''
        for i, res in enumerate(chain):
            if i == 0:
                start_ID = res.id[1]
            try:
                aa = to_one_letter_code[res.resname]
            except KeyError:
                aa= '.'
            sequence += aa
        seqs[chain.id] = sequence
        seqs['%s_st_ID' %chain.id] = start_ID
    return seqs
