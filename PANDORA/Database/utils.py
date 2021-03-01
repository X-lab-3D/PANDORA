#!/usr/bin/python

import os
import csv
import pickle

import numpy as np
from Bio.PDB import PDBParser
from Bio.Data.SCOPData import protein_letters_3to1 as to_one_letter_code

import PANDORA
from PANDORA.Contacts import Contacts

'''
This file contains the following functions:
    - align_peptides(seq1, anch1_seq1, anch2_seq1, seq2, anch1_seq2, anch2_seq2)
    - get_peptides_from_csv(pepts_filename, pept_clmn, allele_clmn, delimiter, skip_first_line = True)
    - get_peptides_w_star_from_csv(pepts_filename, pept_clmn, allele_clmn, delimiter, skip_first_line = True)
    - remove_error_templates(pklfile)
    - change_sep_in_ser(indir)
    - move_uncommon_pdbf(indir, outdir, delete_files = False)
    - small_molecule_MHCI(inputfile, chain_M='M', chain_P='P')
    - get_pdb_seq(IDs)
    - get_seqs(pdbf)
    - get_contacts_matrix(contactfile)
    - get_anchors_pMHCII(Bio.PDB object)
    - pocket_residues_pMHCII(Bio.PDB object)
    - get_anchors_pMHCI(Bio.PDB object)
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
    with open(PANDORA.PANDORA_data +'/csv_pkl_files/%s' %pklfile, 'rb') as inpkl:
        IDD = pickle.load(inpkl)
        bad_IDs = pickle.load(inpkl)
        inpkl.close()
    with open(PANDORA.PANDORA_data +'/csv_pkl_files/error_templates.tsv') as infile: #To be removed later
        r = csv.reader(infile, delimiter='\t')
        for i, row in enumerate(r):
            if i != 0:
                try:
                    del IDD[row[0]]
                    bad_IDs[row[0]] = row[1]
                    os.system('mv '+ PANDORA.PANDORA_data +'/PDBs/pMHCI/'+row[0]+'_MP.pdb ' 
                              + PANDORA.PANDORA_data +'/PDBs/unused_templates/')
                except:
                    pass
        infile.close()
    with open(PANDORA.PANDORA_data +'/csv_pkl_files/' + pklfile, 'wb') as outpkl:
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

        os.system('mv %s %s' %(outfile_name, indir +'/' + pdbfile))

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

def small_molecule_MHCI(inputfile, chain_M='M', chain_P='P', remove_dist_file = True):
    '''
    @author: Rafaella Buzatu
    Takes as input a file (which can be either a pdb file or a contacts list)
    name of chain M and name of chain P

    Args:
        inputfile(str) : Path to either a pdb file or a contacts list
        chain_M(str) : chain ID of MHC-I alpha chain
        chain_P(str) : chain ID of MHC-I bound peptide
        remove_dist_file(bool) : if True, deletes the distance file generated.

    Returns:
        small_mol (str) : Name of the small molecule in the pocket. If no
                          molecules are present, returns None
    #function outputs the name of the small molecule if there is one inside the binding pocket; otherwise, returns None
    '''

    cutoff =str(6)
    pocket = [ 7, 9, 26, 28, 1007, 1009, 1024, 1026]

    ### if the input is a pdb file, the contacts file is calculated
    if inputfile.endswith('.pdb'):
        contactfile = PANDORA.PANDORA_data +'/dist_files/all_contacts_'+ inputfile.split('/')[-1].split('.')[0] +'.list'
        os.popen(PANDORA.PANDORA_path +'/tools/contact-chainID_allAtoms '+ inputfile +' '+ cutoff +' > '+ contactfile).read()

    ### if the contacts file list is used as input, it is used further in the script
    elif inputfile.endswith('.list'):
        contactfile = inputfile

    with open(contactfile) as contacts:

        molecule_in_pocket = []  #list where the name of small molecule(s) is/are appended when found in contact list
        count_molecule = [ 0, 0, 0, 0, 0]  #counter for the frequency of contacts between the small moelcule of corresponding index and the pocket

        ### For cases where numbering starts from 1000, we use the updated pocket residues
        line1=contacts.readline()
        m_aa_id_check = line1.split("\t")[2]
        if int(m_aa_id_check) >1000:
            for i in pocket:
                i = i+1000

        for line in contacts:
            cm_aa_id = line.split("\t")[1]
            m_aa_id = line.split("\t")[2]
            cp_aa_id = line.split("\t")[6]
            molecule = line.split("\t")[5]

            ### Check for contacts that are not between the MHC chain
            ### and peptide or water molecules
            if cm_aa_id == chain_M:
                if cp_aa_id != chain_P  and cp_aa_id != 'B' and molecule != 'HOH':

                    ### Once a suitable contact is found, molecule name is appended in
                    ### moelcule_in_pocket if not previously there
                    ### and the frequency count increases
                    for aa in pocket:
                        try:
                            if int(m_aa_id) == aa:

                                if molecule in molecule_in_pocket:
                                    count_molecule [molecule_in_pocket.index(molecule)] += 1
                                else:
                                    molecule_in_pocket.append(molecule);
                                    count_molecule [molecule_in_pocket.index(molecule)] = 1

                        except ValueError:
                            continue

    if remove_dist_file == True:
        os.system('rm %s ' %contactfile)
        
    if np.max(count_molecule)>0:
        small_mol = molecule_in_pocket[np.argmax(count_molecule)]
        return small_mol
    else:
        return None


def get_pdb_seq(IDs):
    sequences = []
    empty_seqs = []
    for ID in IDs:
        pdbf = PANDORA.PANDORA_data +'/PDBs/pMHCI/%s_MP.pdb' %ID
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

def get_contacts_dict(contactfile):
    temp_contacts_dict = {}
    chains = []
    with open(contactfile) as infile:
        for line in infile:
            donor = line.split('\t')[1]
            acceptor  = line.split('\t')[6]
            if donor not in chains and donor != ' ':
                chains.append(donor)
            if acceptor not in chains and acceptor != ' ':
                chains.append(acceptor)
            if donor != ' ' and acceptor != ' ':
                try: 
                    temp_contacts_dict[(donor, acceptor)] += 1
                    temp_contacts_dict[(acceptor, donor)] += 1
                except KeyError:
                    temp_contacts_dict[(donor, acceptor)] = 1
                    temp_contacts_dict[(acceptor, donor)] = 1
    
    contacts_dict = {}
    for chain in chains:
        contacts_dict[chain] = {}
        for other_chain in chains: 
            if other_chain != chain:
                try:
                    contacts_dict[chain][other_chain] = temp_contacts_dict[chain, other_chain]
                except KeyError:
                    try:
                        contacts_dict[chain][other_chain] = temp_contacts_dict[other_chain, chain]
                    except KeyError:
                        contacts_dict[chain][other_chain] = 0
        #contacts_dict[chain] = {key : 0 for key in chains if key != chain}
        #try:
        #    temp_contacts_dict[chain, ]
    return contacts_dict


def get_anchors_pMHCII(pdb):
    ''' Find the peptide anchor positions of a pMHC2 template structure.

    Args:
        pdb: Bio.PDB object

    Returns: (list) anchor positions of anchor 1,2,3 and 4

    '''

    # todo: This code is very verbose. The anchor finding process can probably be written in a more comprehensive way.
    # uses the function below  pocket_residues () that adapts the pocket residue to the PDB numbering
    pocket_M, pocket_N = pocket_residues_pMHCII(pdb)
    cutoff = 18

    # Calculate contacts
    contacts = Contacts.Contacts(pdb, cutoff=cutoff).chain_contacts

    # vectors that will be used to count the frequency; each value represents the contact frequency between the
    # previously defined pocket residues and the peptide residue corresponding to the vector index
    contact_1 = [0] * 30
    contact_2 = [0] * 30
    contact_3 = [0] * 30
    contact_4 = [0] * 30

    # read each line in the contacts file
    # check for distances between the MHCII chains and the peptide that are smaller than values specified in the pocket_M/N dictionaries
    for line in contacts:
        chain_1 = line[1]
        chain_2 = line[5]
        # print(chain_1, chain_2)

        # the folllowing if statement checks the order of P, M in the contact file
        if chain_1 == 'P':
            cp_aa_id = chain_1
            cm_aa_id = chain_2
            m_aa_id = line[6]
            p_aa_id = line[2]
            atom_name = line[3]
            distance = line[-1]
        else:
            cm_aa_id = chain_1
            cp_aa_id = chain_2
            m_aa_id = line[2]
            p_aa_id = line[6]
            atom_name = line[7]
            distance = line[-1]

        # checks only the lines concerning the C alpha atom of the peptide residues/ C1 atom of modified residues
        if cp_aa_id == 'P':
            if atom_name == 'CA' or atom_name == 'C1':
                # Checks pocket in chain M
                if cm_aa_id == 'M':
                    n = 0
                    for aa in pocket_M.get('pocket_anch1'):
                        try:
                            if int(m_aa_id) == aa:
                                if (float(distance) <= pocket_M.get('pocket_1_distance')[n]):
                                    contact_1[int(p_aa_id)] += 1
                        except ValueError:
                            continue
                        n += 1

                    n = 0
                    for aa in pocket_M.get('pocket_anch3'):
                        try:
                            if int(m_aa_id) == aa:
                                if (float(distance) <= pocket_M.get('pocket_3_distance')[n]):
                                    contact_3[int(p_aa_id)] += 1
                        except ValueError:
                            continue
                        n += 1

                # Check pocket in chain N
                elif cm_aa_id == 'N':
                    n = 0
                    for aa in pocket_N.get('pocket_anch2'):
                        try:
                            if int(m_aa_id) == aa:
                                if (float(distance) <= pocket_N.get('pocket_2_distance')[n]):
                                    contact_2[int(p_aa_id)] += 1

                        except ValueError:
                            continue
                        n += 1

                    n = 0
                    for aa in pocket_N.get('pocket_anch3'):
                        try:
                            if int(m_aa_id) == aa:
                                if (float(distance) <= pocket_N.get('pocket_3_distance')[n]):
                                    contact_3[int(p_aa_id)] += 1

                        except ValueError:
                            continue
                        n += 1

                    n = 0
                    for aa in pocket_N.get('pocket_anch4'):
                        try:
                            if int(m_aa_id) == aa:
                                if (float(distance) <= pocket_N.get('pocket_4_distance')[n]):
                                    contact_4[int(p_aa_id)] += 1
                        except ValueError:
                            continue
                        n += 1

    anchor_1 = np.argmax(contact_1)
    anchor_2 = np.argmax(contact_2)
    anchor_3 = np.argmax(contact_3)
    anchor_4 = np.argmax(contact_4)

    anchors = [anchor_1, anchor_2, anchor_3, anchor_4]
    anchors.sort()

    return anchors


def pocket_residues_pMHCII(pdb):
    '''function takes as input pdbfile and outputs 2 dictionaries that define the peptide binding pocket in chains
    M, N of MHCII

    Args:
        pdb: Bio.PDB object

    Returns:(Tuple of two dicts) of anchor binding positions of MHCII for the M and N chain

    '''

    ###

    # normal pocket  numbering; pocket_anchn defines the residues of chain M/N that form the binding pocket of anchor n
    # pocket_n_distance represents the distances between the nth anchor residue and corresponding nth pocket residues from the first list
    pocket_M = {'pocket_anch1': [28, 33, 31], 'pocket_1_distance': [15, 10.5, 18],
                'pocket_anch3': [10], 'pocket_3_distance': [13.5]}
    pocket_N = {'pocket_anch2': [10, 11, 21, 22, 23], 'pocket_2_distance': [10.7, 10.5, 14.7, 10, 14],
                'pocket_anch3': [6, 7], 'pocket_3_distance': [13, 7.5],
                'pocket_anch4': [28, 27, 31, 33], 'pocket_4_distance': [14, 12.4, 15.7, 13]}

    # the pocket residues are defined differently when the numbering in the M or/and N chain starts from 1000
    pocket_M_renumbered = {'pocket_anch1': [1028, 1033, 1031], 'pocket_1_distance': [15, 10.5, 18],
                           'pocket_anch3': [1010], 'pocket_3_distance': [13.5]}
    pocket_N_renumbered = {'pocket_anch2': [1010, 1011, 1021, 1022, 1023],
                           'pocket_2_distance': [10.7, 10.5, 14.7, 10, 14],
                           'pocket_anch3': [1006, 1007], 'pocket_3_distance': [13, 7.5],
                           'pocket_anch4': [1028, 1027, 1031, 1033], 'pocket_4_distance': [14, 12.5, 15.7, 13]}

    # for each chain (M/N), the first 5 residues are checked; if one of them is >1000,
    # the pocket is defined using the pocket_M/N_renumbered dictionary

    for chain in pdb.get_chains():
        if (chain.get_id() == 'M'):
            n = 0
            for residue in chain.get_residues():
                n += 1
                if n <= 7:
                    try:
                        if int(residue.id[1] > 1000):
                            pocket_M = pocket_M_renumbered
                    except ValueError:
                        continue
                else:
                    break

        if (chain.get_id() == 'N'):
            n = 0
            for residue in chain.get_residues():
                n += 1
                if n <= 7:
                    try:
                        if int(residue.id[1] > 1000):
                            pocket_N = pocket_N_renumbered
                    except ValueError:
                        continue
                else:
                    break

    return (pocket_M, pocket_N)


### ### Function that finds the anchors
def get_anchors_pMHCI(pdb):
    ''' Find the peptide anchor residues of a pMHCI structure using the Contacts class

    Args:
        pdb: Bio.PDB object
    Returns: tuple of anchor residues (int)

    '''

    cutoff = 11.2
    ###Define pocket residues
    pocket = { 'pocket_anch1' : [24, 25, 35, 36],
              'pocket_anch2' : [81, 1005, 1027, 1028, 1029, 1030, 1033]}
    pocket_wrong_numbering = { 'pocket_anch1' : [1024, 1025, 1035, 1036],
                    'pocket_anch2' : [1081, 2005, 2027, 2028, 2029, 2030, 2033]}

    # Template contacts
    contacts = Contacts.Contacts(pdb, cutoff=cutoff).chain_contacts

    contact_1 = [0]*25
    contact_2 = [0]*25

    #for cases where numbering starts from 1000, we use the updated pocket residues
    if int(contacts[0][6]) >1000:
        pocket = pocket_wrong_numbering

    for line in contacts:
        #looks at one residue at a time
        m_aa_id = line[6]
        p_aa_id = line[2]
        atom_name = line[3]

        if atom_name == 'CA':
            for aa in pocket.get('pocket_anch1'):
                try:
                    if int(m_aa_id) == aa:
                        contact_1[int(p_aa_id)]+=1
                except ValueError:
                    continue

            for aa in pocket.get('pocket_anch2'):
                try:
                    if int(m_aa_id) == aa:
                        contact_2[int(p_aa_id)]+=1
                except ValueError:
                    continue

    anchor_1 = np.argmax(contact_1)
    anchor_2 = np.argmax(contact_2)

    return (anchor_1, anchor_2)