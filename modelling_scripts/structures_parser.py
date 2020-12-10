#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 22:54:08 2020

@author: dario
"""
import os
import gzip
import shutil

### TODO: Imports to be changed
from data_prep import get_seqs
from modelling_scripts import del_uncommon_residues_pdbs as durp


def parse_pMHCI_pdbs(ids_list):
    
    ### Preparation
    
    # Variables
    IDs = ids_list
    cwd = os.getcwd()
    indir = 'data/PDBs/IMGT_retrieved/IMGT3DFlatFiles'
    outdir = 'data/PDBs/pMHCI'
    unused_pdbs_dir = 'data/PDBs/unused_templates'
    
    bad_IDs = {}
    err_1 = 0
    err_2 = 0
    err_3 = 0
    err_4 = 0
    err_5 = 0
    err_6 = 0
    err_7 = 0

    # Distance c++ code
    if "contact-chainID_allAtoms" not in os.listdir('%s/modelling_scripts' %cwd):
        os.popen('g++ %s/modelling_scripts/contact-chainID_allAtoms.cpp -o %s/modelling_scripts/contact-chainID_allAtoms' %(cwd, cwd)).read()
    
    
    IDs_dict = {}
    for entry in IDs:
        ID = entry[0]
        ID = ID.upper()
        ID = ID.rstrip()

        with gzip.open('%s/IMGT-%s.pdb.gz' %(indir, ID), 'rb') as f_in:
            with open('%s/%s.pdb' %(outdir, ID), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                
        #Delete .gz file
        #os.popen('rm %s/%s.pdb.gz' %(outdir, ID)).read()

        original_filepath = '%s/%s.pdb' %(outdir, ID)
        selected_chains_filepath = '%s/%s_AC.pdb'%(outdir, ID)
        a_renamed_filepath = '%s/%s_MC.pdb'%(outdir, ID)
        new_filepath = '%s/%s_MP.pdb'%(outdir, ID)
        header_filepath = '%s/%s_header.pdb'%(outdir, ID)
        onlyM_filepath = '%s/%s_M.pdb'%(outdir, ID)
        onlyP_filepath = '%s/%s_P.pdb'%(outdir, ID)

        print('Parsing %s' %ID)
        aID = False
        pID = False
        count_dict = {}

        try:
            seqs = get_seqs(original_filepath)
        except (IndexError, ValueError):
            bad_IDs[ID] = '#1 Parser'
            print('ERROR TYPE #1: Something went wrong with PDBParser. %s added to Bad_IDs' %ID)
            print('##################################')
            err_1 += 1
            os.system('mv %s/%s.pdb %s/parsing_error/ '%(outdir, ID ,unused_pdbs_dir))
            continue

        ### Selecting Alpha chain and peptide chain
        ### TODO: parse REMARK info to get chain IDs and allele information
        
        putative_pID = {}
        for chain in seqs:
            if type(seqs[chain]) == str:
                length = len(seqs[chain])
                if chain != ' ':               #This might be removed, since get_seqs() is  not outputting anymore with empty chain ids
                    count_dict[chain] = length
                if length > 250 and length < 300 and not aID:
                    aID = chain
                elif length > 7 and length < 25:
                    putative_pID[chain] = 0
        if len(putative_pID) == 1:
            pID = list(putative_pID.keys())[0]
        elif len(putative_pID) > 1:
            os.system('./modelling_scripts/contact-chainID_allAtoms %s 5 > ./data/dist_files/%s.dist' %(original_filepath, ID))
            contacts = []
            for line in open('./data/dist_files/%s.dist' %ID, 'r'):
                contact = (line.split('\t')[1], line.split('\t')[6])
                contacts.append(contact)
            for candidate in putative_pID:
                for contact in contacts:
                    if contact[0] == aID:
                        if contact[1] == candidate:
                            putative_pID[candidate] += 1
            pID = (sorted(putative_pID.items(), key=lambda x: x[1], reverse = True))[0][0]

        # Checking chainIDs
        try:
            if aID.islower() or pID.islower():
                bad_IDs[ID] = '#2 Lowercase_chainID'
                print('ERROR TYPE #2: Lower case chain ID. %s added to Bad_IDs' %ID)
                print('##################################')
                err_2 += 1
                os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
                continue
            else:
                pass
        except AttributeError:
            bad_IDs[ID] = '#3'
            print('ERROR TYPE #3: False in chain ID. %s added to Bad_IDs' %ID)
            print(count_dict)
            print('##################################')
            err_3 += 1
            os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
            continue


        ###### Let's classify the structures and understand if they are good or not.
        pept_count = 0
        Amhc = 0
        Bmhc = 0
        Atcr = 0
        Btcr = 0
        for key in count_dict:
            if count_dict[key] > 7 or count_dict[key] < 25:
                pept_count += 1
            elif count_dict[key] > 85 or count_dict[key] < 109:
                Bmhc +=1
            elif count_dict[key] > 190 or count_dict[key] < 210:
                Atcr += 1
            elif count_dict[key] > 209 or count_dict[key] < 250:
                Btcr += 1
            elif count_dict[key] > 250 or count_dict[key] < 300:
                Amhc += 1
            else:
                bad_IDs[ID] = '#4 Unrecognized_chain_lenght'
                print('ERROR TYPE #4: Unrecognized chain lenght. %s added to Bad_IDs' %ID)
                print('##################################')
                err_4 += 1
                os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
                continue
        '''
        if Atcr != 0 or Btcr != 0:
            bad_IDs[ID] = '#5 TCR_inside'
            print('ERROR TYPE #5: TCRs in structure. %s added to Bad_IDs' %ID)
            print('##################################')
            err_5 += 1
            os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
            continue
        
        ######
        else:
        '''
        count_dict['allele'] = entry[1]
        count_dict['pept_seq'] = seqs[pID]
        IDs_dict[ID] = count_dict
        os.system('pdb_selchain -%s,%s %s > %s' %(aID, pID, original_filepath, selected_chains_filepath))
        if pID == 'M':
            os.system('pdb_rplchain -%s:P %s > %s' %(pID, selected_chains_filepath, a_renamed_filepath))
            os.system('pdb_rplchain -%s:M %s > %s' %(aID, a_renamed_filepath, new_filepath))
        else:
            os.system('pdb_rplchain -%s:M %s > %s' %(aID, selected_chains_filepath, a_renamed_filepath))
            os.system('pdb_rplchain -%s:P %s > %s' %(pID, a_renamed_filepath, new_filepath))
        try:
            os.system('rm %s' %original_filepath)
        except:
            pass
        try:
            os.system('rm %s' %selected_chains_filepath)
        except:
            pass
        try:
            os.system('rm %s' %a_renamed_filepath)
        except:
            pass

    ### ADD CORRECTION FOR SEP, F2F, etc.
    os.system('python ./tools/change_sep_in_ser.py %s' %outdir)
    print('Removing uncommon residue files')
    uncommon_pdbs = durp.del_uncommon_pdbf(outdir, unused_pdbs_dir + '/non_canonical_res')
    for u_pdb in uncommon_pdbs:
        try:
            del IDs_dict[u_pdb]
            bad_IDs[u_pdb] = '#6 deleting_error'
            err_6 += 1
        except KeyError:
            print("Tried to delete %s from dataset error #6, error encountered" %u_pdb)

    ### New Errors to be fixed:
    ### #8: Residue 1A
    ###
    ### New Errors to be removed from dataset:
    ### #9: Atipycal anchor position (e.g. P1 instead of P2)
    ### #10: Small Molecule in the binding pocket

    bad_IDs['errors'] = {'#1': err_1, '#2': err_2, '#3': err_3,
                         '#4': err_4, '#5': err_5, '#6': err_6, '#7': err_7}
    
    IDd = open("data/csv_pkl_files/IDs_ChainsCounts_dict.pkl", "wb")
    pickle.dump(IDs_dict, IDd)
    pickle.dump(bad_IDs, IDd)
    IDd.close()
    
    print('BAD IDs:')
    print(bad_IDs)
    '''
    with open('data/csv_pkl_files/IDs_dict.tsv', 'wt') as outfile:
        tsv_writer = csv.writer(outfile, delimiter='\t')
        tsv_writer.writerow(['PDB_ID', 'CHAIN_ID', 'LENGTH', 'ALLELE'])
        for ID in IDs_dict:
            for chain_ID in IDs_dict[ID]:
                if chain_ID != 'allele':
                    tsv_writer.writerow([ID, chain_ID, IDs_dict[ID][chain_ID], IDs_dict[ID]['allele']])
    '''
    #TODO: clean dist_files folder
    remove_error_templates()
    return IDs_dict, bad_IDs
    pass

def parse_pMHCII_pdbs():
    pass

def parse_TCRpMHCI_pdbs():
    pass

def parse_TCRpMHCII_pdbs():
    pass