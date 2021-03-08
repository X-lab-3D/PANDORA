#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 22:54:08 2020

@author: Dario Marzella
"""
import os
import gzip
import shutil
import pickle
from copy import deepcopy
from operator import xor

import PANDORA
import PANDORA.junk.utils as utils
from Bio.PDB import PDBParser


def get_chainid_alleles_MHCI(pdbf):
    '''
    Takes as input an IMGT preprocessed PDB file of p:MHC I.
    Returns a dictionary containing alleles andrelative identity scores for each
    G-domain in the given pdb from the REMARK.

    Args:
        pdbf(str) : path to IMGT pdb file
    '''
    #test: multiple chains 3GJG, multiple alleles 1AO7
    ### Parsing file and extracting remarks
    with open(pdbf) as infile:
        remarks = []
        for line in infile:
            if line.startswith('REMARK 410'):
                row = [x for x in line.rstrip().split(' ') if x != '']
                del row[:2]
                remarks.append(row)
    remarks = [x for x in remarks if x != []]

    ### Dividing each remark section into a chains dictionary
    chains = {}
    flag = False
    for row in remarks:
        if row[0] == 'Chain' and row[1] == 'ID' and len(row)== 4:
            chainID = row[2][-1]
            chains[chainID] = []
            chains[chainID].append(row)
            flag = True
        elif flag == True:
            chains[chainID].append(row)

    ### Extracting MHC I Alpha chains
    mhc_a = {}  # MHC I Alpha
    for chain in chains:
        try:
            if chains[chain][1][3] == 'I-ALPHA':
                mhc_a[chain] = chains[chain]
        except:
            pass

    ### Extracting alleles info
    mhc_a_alleles = {}
    for chain in mhc_a:
        G_dom_alleles = {'G-ALPHA1': [], 'G-ALPHA2': []}
        key = False
        for row in mhc_a[chain]:
            if row[0] == 'G-DOMAIN':
                try:
                    if row[3] == 'description' and row[4] == 'G-ALPHA1':
                        key = 'G-ALPHA1'
                    elif row[3] == 'description' and row[4] == 'G-ALPHA2':
                        key = 'G-ALPHA2'
                    elif key:
                        if row[2] == 'gene' and row[3] == 'and' and row[4] == 'allele':
                            G_dom_alleles[key] += row[5:]
                        else:
                            key = False
                except IndexError:
                    pass
        mhc_a_alleles[chain] = deepcopy(G_dom_alleles)

    mhc_a_alleles_percs = {}
    for chain in mhc_a_alleles:
        mhc_a_alleles_percs[chain] = {}
        for key in mhc_a_alleles[chain]:
            mhc_a_alleles_percs[chain][key] = {}
            ### Allele info are always given with four elements: Gender, Spieces, Allele, Percentage
            for block in range(int(len(mhc_a_alleles[chain][key])/4)):
                allele = mhc_a_alleles[chain][key][2+(4*block)]
                perc = float(mhc_a_alleles[chain][key][3+(4*block)].replace('(', '').replace('%)', '').replace(',',''))
                mhc_a_alleles_percs[chain][key][allele] = perc

    return mhc_a_alleles_percs

def select_alleles_set_MHCI(chain_alleles_percs):
    '''
    Takes as input an dictionary .
    Returns the allele list and the identity scores for each
    G-domain in the given pdb from the REMARK.

    Args:
        chain_alleles_percs(dict) : output dictionary of get_chainid_alleles_MHCI() for only one chain

    Example:
    >>> alleles = get_chainid_alleles_MHCI('PANDORA_files/data/PDBs/pMHCI/1K5N_MP.pdb')
    >>> selected_alleles = select_alleles_set_MHCI(alleles['M'])

    '''


    if xor(chain_alleles_percs['G-ALPHA1'] == {}, chain_alleles_percs['G-ALPHA2'] == {}):
        if chain_alleles_percs['G-ALPHA1'] != {}:
            non_empty = 'G-ALPHA1'
            empty = 'G-ALPHA2'
        elif chain_alleles_percs['G-ALPHA2'] != {}:
            non_empty = 'G-ALPHA2'
            empty = 'G-ALPHA1'
        top_perc = max(list(chain_alleles_percs[non_empty].values()))
        to_delete = []
        for allele in chain_alleles_percs[non_empty]:
            # Checking if the percentace is the highest
            if chain_alleles_percs[non_empty][allele] < top_perc:
                to_delete.append(allele)
        for allele in to_delete:
            chain_alleles_percs[non_empty].pop(allele)
        chain_alleles_percs.pop(empty)

    elif chain_alleles_percs['G-ALPHA1'] == {} and chain_alleles_percs['G-ALPHA2'] == {}:
        return None
    else:
        for i, domain in enumerate(chain_alleles_percs):
            other_domain = list(chain_alleles_percs.keys())[not i]

            top_perc = max(list(chain_alleles_percs[domain].values()))
            to_delete = []
            for allele in chain_alleles_percs[domain]:
                # Checking if the percentace is the highest
                if chain_alleles_percs[domain][allele] < top_perc:
                    to_delete.append(allele)
                # Checking if the allele is present in both domains
                if allele not in chain_alleles_percs[other_domain]:
                    to_delete.append(allele)
            for allele in to_delete:
                chain_alleles_percs[domain].pop(allele)

    selected_alleles = [chain_alleles_percs[dom] for dom in chain_alleles_percs]
    selected_alleles = [list(x.keys()) for x in selected_alleles]
    selected_alleles = list(set(sum(selected_alleles, [])))

    return selected_alleles


def parse_pMHCI_pdbs(ids_list, out_pkl = 'IDs_and_bad_IDs_dict.pkl', remove_dist_file = True):
    '''
    Parses pMHCI structures for PANDORA in the following steps for each structure:
    1. Uncompress and move to pMHCI folder (outdir).
    2. Indetify one Alpha Chain and its Peptide
    3. Move Alpha Chain and Peptide to a new PDB file, name <ID>_MP.pdb containing only chain M (alpha chain) and P (peptide)
    4. Modify all phosphorilates serines in serines
    5. Remove from the clean dataset every structure containing non-canonical residues

    '''

    ### Preparation

    # Variables
    IDs = ids_list
    indir = PANDORA.PANDORA_data + '/PDBs/IMGT_retrieved/IMGT3DFlatFiles'
    outdir = PANDORA.PANDORA_data + '/PDBs/pMHCI'
    unused_pdbs_dir = PANDORA.PANDORA_data + '/PDBs/unused_templates'

    bad_IDs = {}
    err_1 = 0
    err_2 = 0
    err_3 = 0
    err_4 = 0
    err_5 = 0
    err_6 = 0
    err_7 = 0

    IDs_dict = {}
    for ID in IDs:
        ID = ID.upper()
        ID = ID.rstrip()

        ### Try to unzip the file
        try:
            with gzip.open('%s/IMGT-%s.pdb.gz' %(indir, ID), 'rb') as f_in:
                with open('%s/%s.pdb' %(outdir, ID), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
        except FileNotFoundError:
            bad_IDs[ID] = '#1 FileNotFound'
            print('ERROR TYPE #1: File not found. %s added to Bad_IDs' %ID)
            print('##################################')
            err_1 += 1
            os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
            continue

        ### Paths to use later
        original_filepath = '%s/%s.pdb' %(outdir, ID)
        selected_chains_filepath = '%s/%s_AC.pdb'%(outdir, ID)
        a_renamed_filepath = '%s/%s_MC.pdb'%(outdir, ID)
        new_filepath = '%s/%s_MP.pdb'%(outdir, ID)
        header_filepath = '%s/%s_header.pdb'%(outdir, ID)
        onlyM_filepath = '%s/%s_M.pdb'%(outdir, ID)
        onlyP_filepath = '%s/%s_P.pdb'%(outdir, ID)
        persian_filepath = '%s/%s_persian.pdb'%(outdir, ID)

        print('Parsing %s' %ID)
        aID = False
        pID = False
        count_dict = {}

        ### Get the sequences
        try:
            seqs = utils.get_seqs(original_filepath)
        except (IndexError, ValueError):
            bad_IDs[ID] = '#2 Parser'
            print('ERROR TYPE #2: Something went wrong with PDBParser. %s added to Bad_IDs' %ID)
            print('##################################')
            err_2 += 1
            os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
            continue

        ### Select Alpha chain and peptide chain

        alpha_chains = get_chainid_alleles_MHCI(original_filepath)

        ### Identify putative peptide and alpha chain
        ### TODO: update this step to make it more flexible
        putative_pID = {}
        for chain in seqs:
            if type(seqs[chain]) == str:
                length = len(seqs[chain])
                count_dict[chain] = length
                if (chain in list(alpha_chains.keys())): # and (not aID) and (length > 250 and length < 300)
                    aID = chain                          # TODO: Try another putative aID if the first does not work
                elif length > 7 and length < 25:
                    putative_pID[chain] = 0

        ### Check which putative peptide is close enough to the alpha chain
        if len(putative_pID) == 1:
            pID = list(putative_pID.keys())[0]
        elif len(putative_pID) > 1:
            dist_file = PANDORA.PANDORA_data + '/dist_files/'+ ID +'.list'
            os.system(PANDORA.PANDORA_path + '/tools/contact-chainID_allAtoms ' + original_filepath + ' 5 > ' + 
                      dist_file)
            contacts = []
            for line in open(dist_file, 'r'):
                contact = (line.split('\t')[1], line.split('\t')[6])
                contacts.append(contact)
            for candidate in putative_pID:
                for contact in contacts:
                    if contact[0] == aID:
                        if contact[1] == candidate:
                            putative_pID[candidate] += 1
            pID = (sorted(putative_pID.items(), key=lambda x: x[1], reverse = True))[0][0]

            if remove_dist_file == True:
                os.system('rm %s/dist_files/%s.list' %(PANDORA.PANDORA_data, ID))
                
        ### Check chainIDs (lowercase or empty chainIDs would raise errors later)
        try:
            if aID.islower() or pID.islower():
                bad_IDs[ID] = '#3 Lowercase_chainID'
                print('ERROR TYPE #3: Lower case chain ID. %s added to Bad_IDs' %ID)
                print('##################################')
                err_3 += 1
                os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
                continue
            else:
                pass
        except AttributeError:
            bad_IDs[ID] = '#4 chainID'
            print('ERROR TYPE #4: False in chain ID. %s added to Bad_IDs' %ID)
            print(count_dict)
            print('##################################')
            err_4 += 1
            os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
            continue


        ### Check sequences lengths to identify molecules inside (obsolete at the moment)
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


        ### Retrive alleles
        alleles = select_alleles_set_MHCI(alpha_chains[aID])
        if alleles == None:
            bad_IDs[ID] = '#5 No_alleles'
            print('ERROR TYPE #5: No available alleles. %s added to Bad_IDs' %ID)
            print('##################################')
            err_5 += 1
            os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
            continue
        else:
            pass


        count_dict['allele'] = alleles
        count_dict['pept_seq'] = seqs[pID]
        IDs_dict[ID] = count_dict

        ### Rename M and P chain, using an intermediate with persian M and P (in order to avoid ID overlaps)
        os.system('pdb_rplchain -%s:پ %s > %s' %(pID, original_filepath, a_renamed_filepath))
        os.system('pdb_rplchain -%s:م %s > %s' %(aID, a_renamed_filepath, selected_chains_filepath))

        os.system('pdb_selchain -م,پ %s > %s' %(selected_chains_filepath, persian_filepath))

        os.system('pdb_rplchain -پ:P %s > %s' %(persian_filepath, a_renamed_filepath))
        os.system('pdb_rplchain -م:M %s > %s' %(a_renamed_filepath, new_filepath))

        ### Rorder the PDB file if chain P comes before chain M
        if list(seqs.keys()).index(pID) < list(seqs.keys()).index(aID):

            ### Write only-header PDB file
            with open(new_filepath, 'r') as inpdb:
                with open(header_filepath, 'w') as outheader:
                    for line in inpdb:
                        if line.startswith('ATOM'):
                            break
                        else:
                            outheader.write(line)
                    outheader.close()
                inpdb.close()

            ### Extract M and P chains in two different PDBs
            os.system('pdb_selchain -M %s > %s' %(new_filepath, onlyM_filepath))
            os.system('pdb_selchain -P %s > %s' %(new_filepath, onlyP_filepath))

            ### Merge the three PDB files together
            with open(new_filepath, 'w') as outpdb:
                with open(header_filepath, 'r') as inheader:
                    for line in inheader:
                        outpdb.write(line)
                    inheader.close()
                with open(onlyM_filepath, 'r') as inMpdb:
                    for line in inMpdb:
                        if line.startswith('ATOM') or line.startswith('TER'):
                            outpdb.write(line)
                    inMpdb.close()
                with open(onlyP_filepath, 'r') as inPpdb:
                    for line in inPpdb:
                        if line.startswith('ATOM') or line.startswith('TER'):
                            outpdb.write(line)
                    inPpdb.close()
                outpdb.write('END')
                outpdb.close()

            ### Remove temporary files
            try:
                os.system('rm %s' %header_filepath)
            except:
                pass
            try:
                os.system('rm %s' %onlyM_filepath)
            except:
                pass
            try:
                os.system('rm %s' %onlyP_filepath)
            except:
                pass

        ### Remove temporary files
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
        try:
            os.system('rm %s' %persian_filepath)
        except:
            pass

        ### Check if there are small molecules in the binding cleft
        try:
            small_molecule = utils.small_molecule_MHCI(new_filepath)
        except IndexError:
            del IDs_dict[ID]
            bad_IDs[ID] = '#4 chainID'
            print('ERROR TYPE #4: False in chain ID. %s added to Bad_IDs' %ID)
            print(count_dict)
            print('##################################')
            err_4 += 1
            os.system('mv %s %s/parsing_errors/ '%(new_filepath,unused_pdbs_dir))
            continue

        if small_molecule != None:
            del IDs_dict[ID]
            bad_IDs[ID] = '#6 Small Molecule'
            print('ERROR TYPE #6: Small Molecule in the binding cleft. %s added to Bad_IDs' %ID)
            print('##################################')
            err_6 += 1
            os.system('mv %s %s/parsing_errors/ '%(new_filepath ,unused_pdbs_dir))
            continue

    ### Correct SEP, F2F, etc. residues into normal residues.
    utils.change_sep_in_ser(outdir)

    print('Removing uncommon residue files')
    uncommon_pdbs = utils.move_uncommon_pdbf(outdir + '/', unused_pdbs_dir + '/non_canonical_res')
    for u_pdb in uncommon_pdbs:
        try:
            del IDs_dict[u_pdb]
            bad_IDs[u_pdb] = '#7 uncommon residue'
            print('ERROR TYPE #7: Uncommon residue PDB. %s added to Bad_IDs' %ID)
            print('##################################')
            err_7 += 1
        except KeyError:
            print("Tried to delete %s from dataset error #7, error encountered" %u_pdb)

    ### New Errors to be fixed:
    ### #8: Residue 1A

    bad_IDs['errors'] = {'#1': err_1, '#2': err_2, '#3': err_3,
                         '#4': err_4, '#5': err_5, '#6': err_6, '#7': err_7}

    if out_pkl:
        IDd = open(PANDORA.PANDORA_data + "/csv_pkl_files/%s" %out_pkl, "wb")
        pickle.dump(IDs_dict, IDd)
        pickle.dump(bad_IDs, IDd)
        IDd.close()

    print('BAD IDs:')
    print(bad_IDs)

    #TODO: clean dist_files folder

    if out_pkl:
        utils.remove_error_templates(out_pkl)

    return IDs_dict, bad_IDs

def get_chainid_alleles_MHCII(pdbf):
    '''
    Takes as input an IMGT preprocessed PDB file of p:MHC II.
    Returns a dictionary containing alleles andrelative identity scores for each
    G-domain in the given pdb from the REMARK.

    Args:
        pdbf(str) : path to IMGT pdb file
    '''
    #test: multiple chains 3GJG, multiple alleles 1AO7
    ### Parsing file and extracting remarks
    with open(pdbf) as infile:
        remarks = []
        for line in infile:
            if line.startswith('REMARK 410'):
                row = [x for x in line.rstrip().split(' ') if x != '']
                del row[:2]
                remarks.append(row)
    remarks = [x for x in remarks if x != []]

    ### Dividing each remark section into a chains dictionary
    chains = {}
    flag = False
    for row in remarks:
        if row[0] == 'Chain' and row[1] == 'ID' and len(row)== 4:
            chainID = row[2][-1]
            chains[chainID] = []
            chains[chainID].append(row)
            flag = True
        elif flag == True:
            chains[chainID].append(row)

    ### Extracting MHC II Alpha and Beta chains
    mhc_a = {}  # MHC II Alpha
    mhc_b = {}  # MHC II Beta
    for chain in chains:
        try:
            if chains[chain][1][3] == 'II-ALPHA':
                mhc_a[chain] = chains[chain]
            elif chains[chain][1][3] == 'II-BETA':# or chains[chain][1][3] == 'II-BETA':
                mhc_b[chain] = chains[chain]
        except:
            pass

    ### Extracting alleles info
    mhc_a_alleles = {x:[] for x in mhc_a}
    mhc_b_alleles = {x:[] for x in mhc_b}
    
    for chain in mhc_a:
        key = False
        for row in mhc_a[chain]:
            if row[0] == 'G-DOMAIN':
                try:
                    if row[3] == 'description' and row[4] == 'G-ALPHA':
                        key = 'G-ALPHA'
                    elif key:
                        if row[2] == 'gene' and row[3] == 'and' and row[4] == 'allele':
                            mhc_a_alleles[chain] += row[5:]
                        else:
                            key = False
                except IndexError:
                    pass

    for chain in mhc_b:
        key = False
        for row in mhc_b[chain]:
            if row[0] == 'G-DOMAIN':
                try:
                    if row[3] == 'description' and row[4] == 'G-BETA':
                        key = 'G-BETA'
                    elif key:
                        if row[2] == 'gene' and row[3] == 'and' and row[4] == 'allele':
                            mhc_b_alleles[chain] += row[5:]
                        else:
                            key = False
                except IndexError:
                    pass

    '''
    mhc_a_alleles_percs = {}
    for chain in mhc_a_alleles:
        mhc_a_alleles_percs[chain] = {}
        for key in mhc_a_alleles[chain]:
            mhc_a_alleles_percs[chain][key] = {}
            ### Allele info are always given with four elements: Gender, Spieces, Allele, Percentage
            for block in range(int(len(mhc_a_alleles[chain][key])/4)):
                allele = mhc_a_alleles[chain][key][2+(4*block)]
                perc = float(mhc_a_alleles[chain][key][3+(4*block)].replace('(', '').replace('%)', '').replace(',',''))
                mhc_a_alleles_percs[chain][key][allele] = perc
    '''
    return {'Alpha': mhc_a_alleles, 'Beta': mhc_b_alleles}

# def parse_pMHCII_pdbs(ids_list, out_pkl = 'pMHC2_IDs_and_bad_IDs_dict.pkl'):
#     '''
#     Parses pMHCII structures for PANDORA in the following steps:
#     1. Uncompress and move every strcuture.
#
#     FOLLOWING STILL TO DEVELOP
#     2. Indetify one Alpha Chain and its Peptide
#     3. Move Alpha Chain and Peptide to a new PDB file, name <ID>_MP.pdb containing only chain M (alpha chain) and P (peptide)
#     4. Modify all phosphorilates serines in serines
#     5. Remove from the clean dataset every structure containing non-canonical residues
#
#     '''
#
#     ### Preparation
#
#     # Variables
#     IDs = ids_list
#     cwd = os.getcwd()
#     indir =  PANDORA.PANDORA_data + '/PDBs/IMGT_retrieved/IMGT3DFlatFiles'
#     outdir =  PANDORA.PANDORA_data + '/PDBs/pMHCII'
#     unused_pdbs_dir =  PANDORA.PANDORA_data + '/unused_templates'
#
#     bad_IDs = {}
#     err_1 = 0
#     err_2 = 0
#     err_3 = 0
#     err_4 = 0
#     err_5 = 0
#     err_6 = 0
#     err_7 = 0
#
#     IDs_dict = {}
#     for ID in IDs:
#         ID = ID.upper()
#         ID = ID.rstrip()
#
#         try:
#             with gzip.open('%s/IMGT-%s.pdb.gz' %(indir, ID), 'rb') as f_in:
#                 with open('%s/%s.pdb' %(outdir, ID), 'wb') as f_out:
#                     shutil.copyfileobj(f_in, f_out)
#         except FileNotFoundError:
#             bad_IDs[ID] = '#1 FileNotFound'
#             print('ERROR TYPE #1: File not found. %s added to Bad_IDs' %ID)
#             print('##################################')
#             err_1 += 1
#             os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
#             continue
#
#         original_filepath = '%s/%s.pdb' %(outdir, ID)
#         '''
#         selected_chains_filepath = '%s/%s_AC.pdb'%(outdir, ID)
#         a_renamed_filepath = '%s/%s_MC.pdb'%(outdir, ID)
#         new_filepath = '%s/%s_MP.pdb'%(outdir, ID)
#         header_filepath = '%s/%s_header.pdb'%(outdir, ID)
#         onlyM_filepath = '%s/%s_M.pdb'%(outdir, ID)
#         onlyP_filepath = '%s/%s_P.pdb'%(outdir, ID)
#         persian_filepath = '%s/%s_persian.pdb'%(outdir, ID)
#         '''
#
#         print('Parsing %s' %ID)
#         aID = False
#         bID = False
#         pID = False
#         count_dict = {}
#
#         try:
#             seqs = utils.get_seqs(original_filepath)
#         except (IndexError, ValueError):
#             bad_IDs[ID] = '#2 Parser'
#             print('ERROR TYPE #2: Something went wrong with PDBParser. %s added to Bad_IDs' %ID)
#             print('##################################')
#             err_2 += 1
#             os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
#             continue
#
#         ### Get Alpha and Beta chains info
#         chains = get_chainid_alleles_MHCII(original_filepath)
#         IDs_dict[ID] = chains
#         a_chains = list(chains['Alpha'].keys())
#         b_chains = list(chains['Beta'].keys())
#
#         ### Identify putative peptides
#         putative_pID = []
#         for chain in seqs:
#             if type(seqs[chain]) == str:
#                 length = len(seqs[chain])
#                 if length > 7 and length < 25:
#                     putative_pID.append(chain)
#
#         if len(putative_pID) == 1:
#             pID = list(putative_pID.keys())[0]
#         elif len(putative_pID) > 1:
#             os.system( PANDORA.PANDORA_path + '/tools/contact-chainID_allAtoms %s 5 > ./PANDORA_files/data/dist_files/%s.list' %(original_filepath, ID))
#             contactfile =  PANDORA.PANDORA_data + '/dist_files/%s.list' %ID
#             contacts_dict = utils.get_contacts_dict(contactfile)
#             for a_chain in list(chains['Alpha'].keys()):
#                 ab_contacts = {key : contacts_dict[key] for key in list(contacts_dict.keys())if key in b_chains}
#                 b_chain = max(ab_contacts, key=ab_contacts.get)
#                 if ab_contacts[b_chains] > 200:  # Check if alpha and beta chain have enough contacts.
#                     break
#             # TODO: do the same thing with putative_pID
#
#
#
#         '''
#         putative_pID_contacts = {}
#         ### Identify putative peptide and alpha chain
#         for chain in seqs:
#             if type(seqs[chain]) == str:
#                 length = len(seqs[chain])
#                 count_dict[chain] = length
#                 if (chain in list(chains['Alpha'].keys())): # and (not aID) and (length > 250 and length < 300)
#                     aID = chain                          # TODO: Try another putative aID if the first does not work
#                 elif if (chain in list(chains['Beta'].keys())):
#                     bID = chain
#                 elif length > 7 and length < 25:
#                     putative_pID_contacts[chain] = 0
#
#         ### Check which putative peptide is close enough to the alpha chain
#         if len(putative_pID) == 1:
#             pID = list(putative_pID.keys())[0]
#         elif len(putative_pID) > 1:
#             os.system('./PANDORA/tools/contact-chainID_allAtoms %s 5 > ./PANDORA_files/data/dist_files/%s.dist' %(original_filepath, ID))
#             contacts = []
#             for line in open('./PANDORA_files/data/dist_files/%s.dist' %ID, 'r'):
#                 contact = (line.split('\t')[1], line.split('\t')[6])
#                 contacts.append(contact)
#             for candidate in putative_pID:
#                 for contact in contacts:
#                     if contact[0] == aID:
#                         if contact[1] == candidate:
#                             putative_pID[candidate] += 1
#             pID = (sorted(putative_pID.items(), key=lambda x: x[1], reverse = True))[0][0]
#             '''
#
#     with open( PANDORA.PANDORA_data + '/csv_pkl_files/%s' %out_pkl, 'wb') as outfile:
#         pickle.dump(IDs_dict, outfile)
#         pickle.dump(bad_IDs, outfile)
#
#     return IDs_dict, bad_IDs


def unzip_pdb(ID, indir, outdir):
    """ Unzips a pdb, move it to another directory and return the filepath

    :param indir: location of pdb.gz files
    :param outdir: output location
    :return: filepath of .pdb
    """
    ## unzip pdb and move to outdir
    try:
        with gzip.open('%s/IMGT-%s.pdb.gz' % (indir, ID), 'rb') as f_in:
            # Check if the file is empty
            if f_in.seek(0, whence=2) == 0:
                raise Exception('File is empty')

        with gzip.open('%s/IMGT-%s.pdb.gz' % (indir, ID), 'rb') as f_in:
            with open('%s/%s.pdb' % (outdir, ID), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    except FileNotFoundError:
        print('ERROR TYPE #1: File not found. %s' % ID)

    return '%s/%s.pdb' % (outdir, ID)

def find_chains_MHCII(ID, pdb_file,pept_dist = 10):
    ''' Find the Alpha, Beta and Peptide chain in a pdb file. Also checks if there are heteroatoms in the binding
    cleft and checks for the proper size of the peptide chain.

    :param ID: (string) pdb ID
    :param pdb_file: (string) path to pdb file
    :param pept_dist: (int) distance in angstrom between the peptide chain and putative MHC chains to select them.
    :return: [Alpha_chain, Beta_chain, Peptide_chain] and the pdb object
    '''
    parser = PDBParser(QUIET=True)  # Create a parser object, used to read pdb files
    pdb = parser.get_structure(ID, pdb_file)

    # Remove weird chains
    for i in pdb.get_chains():
        for model in pdb:
            for chain in model:
                if chain.id in [' ', '', '1', '2', '3', '4', '5', '6', '7', '8', 'S']:
                    model.detach_child(chain.id)

    try:
    # Find most likely peptide chain: first chain to be 7 < len(chain) < 25
        pept_chain = []
        for chain in pdb.get_chains():
            if len(chain) > 6 and len(chain) < 26: # Is this chain between 7 and 25?
                heteroatoms = False
                for res in chain:
                    if res.id[0] != " " and res.id[0] != 'W': # Check if a res in this chain is a heteroatom
                        heteroatoms = True
                        break
                if heteroatoms == False: # If all residues are oke, add this to the list of peptide chains
                    pept_chain.append(chain.id)
        pept_chain = pept_chain[0] #Take the first pept chain. If there are multiple, they are probably duplicates

        if len(pept_chain) == 0:
            print('Could not find peptide chain')
            raise Exception
    except:
        print('Could not find peptide chain')
        raise Exception

    # Find peptide-structure contacts
    try:
        chain_cont = []
        for pept_res in pdb[0][pept_chain].get_residues(): #loop through all alpha carbons of the peptide residues
            for chain in pdb.get_chains(): #go through all chains and find residues of that chain that contact the peptide
                if chain.id != pept_chain: #dont measure distance between peptide-peptide residues of course
                    for res in chain:
                        if res.id[0] == ' ': #check if the residue is an amino acid and no heteroatom or water
                            try:
                                if res['CA'] - pept_res['CA'] < pept_dist: #distance between peptide and chain residue alpha carbon
                                    chain_cont.append(chain.id)
                            except:
                                pass
    except:
        print('Could not find peptide chain')
        raise Exception

    # Find the two chains that have the most contacts with the peptide. This should also filter out TCR chains
    if len(set(chain_cont)) >= 2:
        MHC_chains = sorted([ss for ss in set(chain_cont)], key=chain_cont.count, reverse=True)[:2]
        MHC_chains = sorted(MHC_chains) # sort alphabetically
        MHC_chains = MHC_chains + [pept_chain]
    else:
        print('Found >2 MHC chains')
        raise Exception

    for i in pdb.get_chains():
        for model in pdb:
            for chain in model:
                if chain.id not in MHC_chains:
                    model.detach_child(chain.id)

    return MHC_chains, pdb

def replace_chain_names(chains, pdb, replacement_chains=['M', 'N', 'P']):
    ''' Replace chain names by another chain name in a bio.pdb object

    :param chains: (list) chains to replace
    :param pdb: bio.pdb object
    :param replacement_chains: (list) replacement names (in order of chains to replace)
    :return: bio.pdb object with changed chain names
    '''

    for i in chains:
        for chain in pdb.get_chains():
            if chain.id == i:
                chain.id = replacement_chains[chains.index(i)]
    return pdb


def write_pdb(pdb, out_path, get_header_from=False):
    ''' Write bio.pdb object to file, can use the header of the original pdb (bio.pdb cant remember file headers)

    :param pdb: bio.pdb object
    :param out_path: (string)
    :param get_header_from: (string) get the header from another pdb file
    '''

    def get_head_and_remarks(pdb_file):
        ''' Get the head and remarks of an IMGT pdb file

        :param pdb_file: path to pdb file
        :return: (list) list of lines
        '''
        # Count until where the header and remarks last
        x = 0
        last_header_line = 0
        with open(pdb_file) as infile:
            header = []
            for line in infile:
                x += 1
                header.append(line[:-1])
                if line.startswith('REMARK'):
                    if last_header_line == 0:
                        last_header_line = x
                    last_header_line += 1
        header = [x for x in header[:last_header_line - 1] if
                  x != []]  # remove everything after the remarks and whitelines
        return header

    def line_prepender(filename, line):
        ''' Add a line in front of a file

        :param filename: (string) filepath
        :param line: (string) line to prepend
        '''
        with open(filename, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(line.rstrip('\r\n') + '\n' + content)

    # If the original pdb file path is given, use that header and paste it before the ATOM lines
    if get_header_from:
        header = get_head_and_remarks(get_header_from)  # get the header

    # Write pdb
    io = PDBIO()
    io.set_structure(pdb)
    io.save(out_path)

    if get_header_from: # Write the original header to pdb file
        for line in header[::-1]:
            line_prepender(out_path, line)

def parse_pMHCII_pdbs(ids_list):
    ''' Clean all MHCII pdb files that have been downloaded from IMGT

    :param ids_list: (list) list of MHCII PDB IDs. core.Database.IDs_list_MHCII
    :return: Writes all cleaned PDBs to the /PDBs/pMHCII/ dir
    '''
    # keep track of the failed PDBs
    bad_ids = []
    # set paths for in and out directories
    indir = PANDORA.PANDORA_data + '/PDBs/IMGT_retrieved/IMGT3DFlatFiles'
    outdir = PANDORA.PANDORA_data + '/PDBs/pMHCII'

    for ID in ids_list:
        ID = '1A6A'
        print('Parsing %s' %ID)
        try:
            # Unzip file (also check if the file is not empty) and save the path of this file
            pdb_file = unzip_pdb(ID, indir, outdir)
            # Find out which chains are the Alpha, Beta and Peptide chain
            MHC_chains, pdb = find_chains_MHCII(ID, pdb_file)
            # Rename the Alpha, Beta and Peptide chain to M,N,P respectively
            pdb = replace_chain_names(MHC_chains, pdb, ['M', 'N', 'P'])
            # Finally, write the cleaned pdb to the output dir. Keep the header of the original file.
            write_pdb(pdb, '%s/%s.pdb' %(outdir, ID), pdb_file)

        except: # If something goes wrong, append the ID to the bad_ids list
            bad_ids.append(ID)
            pass


# for i in pdb.get_chains:
#     print(i)