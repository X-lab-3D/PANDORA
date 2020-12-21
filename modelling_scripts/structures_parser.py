#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 22:54:08 2020

@author: dario
"""
import os
import gzip
import shutil
import pickle
from copy import deepcopy
from operator import xor

### TODO: Imports to be changed
from modelling_scripts.utils import get_seqs
from modelling_scripts.utils import remove_error_templates
import modelling_scripts.del_uncommon_residues_pdbs as durp

def get_chainid_alleles_MHCI(pdbf):
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
    # test = mhc_a_alleles_percs['A']
    if xor(chain_alleles_percs['G-ALPHA1'] == {}, chain_alleles_percs['G-ALPHA2'] == {}):
        if chain_alleles_percs['G-ALPHA1'] != {}:
            non_empty = 'G-ALPHA1'
            empty = 'G-ALPHA2'
        elif chain_alleles_percs['G-ALPHA2'] != {}:
            non_empty = 'G-ALPHA2'
            empty = 'G-APLHA1'
        top_perc = max(list(chain_alleles_percs[non_empty].values()))
        to_delete = []
        for allele in chain_alleles_percs[non_empty]:
            # Checking if the percentace is the highest
            if chain_alleles_percs[non_empty][allele] < top_perc:
                to_delete.append(allele)
        for allele in to_delete:
            chain_alleles_percs[non_empty].pop(allele)
        chain_alleles_percs.pop(empty)
        
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
        

def parse_pMHCI_pdbs(ids_list, out_pkl = 'IDs_and_bad_IDs_dict.pkl'):
    
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
    for ID in IDs:
        ID = ID.upper()
        ID = ID.rstrip()
        
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
            bad_IDs[ID] = '#2 Parser'
            print('ERROR TYPE #2: Something went wrong with PDBParser. %s added to Bad_IDs' %ID)
            print('##################################')
            err_2 += 1
            os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
            continue

        ### Selecting Alpha chain and peptide chain
        ### TODO: parse REMARK info to get chain IDs and allele information
        
        alpha_chains = get_chainid_alleles_MHCI(original_filepath)
        
        putative_pID = {}
        for chain in seqs:
            if type(seqs[chain]) == str:
                length = len(seqs[chain])
                count_dict[chain] = length
                if (chain in list(alpha_chains.keys())): # and (not aID) and (length > 250 and length < 300)
                    aID = chain                          # TODO: Try another putative aID if the first does not work
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
                bad_IDs[ID] = '#3 Lowercase_chainID'
                print('ERROR TYPE #3: Lower case chain ID. %s added to Bad_IDs' %ID)
                print('##################################')
                err_3 += 1
                os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
                continue
            else:
                pass
        except AttributeError:
            bad_IDs[ID] = '#4'
            print('ERROR TYPE #4: False in chain ID. %s added to Bad_IDs' %ID)
            print(count_dict)
            print('##################################')
            err_4 += 1
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
            '''
            else:
                bad_IDs[ID] = '#4 Unrecognized_chain_lenght'
                print('ERROR TYPE #4: Unrecognized chain lenght. %s added to Bad_IDs' %ID)
                print('##################################')
                err_4 += 1
                os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
                continue
            '''
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
        
        #print(alpha_chains[aID])
        if any([1 for key in alpha_chains[aID] if alpha_chains[aID][key] == {}]):
            bad_IDs[ID] = '#5 No_alleles'
            print('ERROR TYPE #5: No available alleles. %s added to Bad_IDs' %ID)
            print('##################################')
            err_5 += 1
            os.system('mv %s/%s.pdb %s/parsing_errors/ '%(outdir, ID ,unused_pdbs_dir))
            continue
        else:
            alleles = select_alleles_set_MHCI(alpha_chains[aID])
        
        
        count_dict['allele'] = alleles
        count_dict['pept_seq'] = seqs[pID]
        IDs_dict[ID] = count_dict
        os.system('pdb_selchain -%s,%s %s > %s' %(aID, pID, original_filepath, selected_chains_filepath))
        if pID == 'M':
            os.system('pdb_rplchain -%s:P %s > %s' %(pID, selected_chains_filepath, a_renamed_filepath))
            os.system('pdb_rplchain -%s:M %s > %s' %(aID, a_renamed_filepath, new_filepath))
        else:
            os.system('pdb_rplchain -%s:M %s > %s' %(aID, selected_chains_filepath, a_renamed_filepath))
            os.system('pdb_rplchain -%s:P %s > %s' %(pID, a_renamed_filepath, new_filepath))
            
        #IF chain P comes before chain M
        if list(seqs.keys()).index(pID) < list(seqs.keys()).index(aID):
            
            with open(new_filepath, 'r') as inpdb:
                with open(header_filepath, 'w') as outheader:
                    for line in inpdb:
                        if line.startswith('ATOM'):
                            break
                        else:
                            outheader.write(line)
                    outheader.close()
                inpdb.close()
            
            os.system('pdb_selchain -M %s > %s' %(new_filepath, onlyM_filepath))
            os.system('pdb_selchain -P %s > %s' %(new_filepath, onlyP_filepath))
            
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
    os.system('python ./tools/change_sep_in_ser.py %s/' %outdir)
    print('Removing uncommon residue files')
    uncommon_pdbs = durp.move_uncommon_pdbf(outdir + '/', unused_pdbs_dir + '/non_canonical_res')
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
    
    IDd = open("data/csv_pkl_files/%s" %out_pkl, "wb")
    pickle.dump(IDs_dict, IDd)
    pickle.dump(bad_IDs, IDd)

    IDd.close()
    
    print('BAD IDs:')
    print(bad_IDs)
    
    #TODO: clean dist_files folder
    
    #remove_error_templates()
    
    return IDs_dict, bad_IDs

def get_chainid_alleles_MHCII(pdbf):
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
    mhc_a_alleles = {}
    mhc_b_alleles = {}
    for chain in mhc_a:
        key = False
        for row in mhc_a[chain]:
            if row[0] == 'G-DOMAIN':
                try:
                    if row[3] == 'description' and row[4] == 'G-ALPHA':
                        key = 'G-ALPHA'
                    elif key:
                        if row[2] == 'gene' and row[3] == 'and' and row[4] == 'allele':
                            mhc_a_alleles[key] += row[5:]
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
                        key = 'G-ALPHA'
                    elif key:
                        if row[2] == 'gene' and row[3] == 'and' and row[4] == 'allele':
                            mhc_b_alleles[key] += row[5:]
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
    return {'Aplha': mhc_a_alleles, 'Beta': mhc_b_alleles}

def parse_pMHCII_pdbs(ids_list):
    
    ### Preparation
    
    # Variables
    IDs = ids_list
    cwd = os.getcwd()
    indir = 'data/PDBs/IMGT_retrieved/IMGT3DFlatFiles'
    outdir = 'data/PDBs/pMHCII'
    unused_pdbs_dir = 'data/PDBs/unused_templates'
    
    bad_IDs = {}
    err_1 = 0
    
    IDs_dict = {}
    for ID in IDs:
        ID = ID.upper()
        ID = ID.rstrip()
        
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
        
        original_filepath = '%s/%s.pdb' %(outdir, ID)
        
        IDs_dict[ID] = get_chainid_alleles_MHCII(ID)
    
    return IDs_dict, bad_IDs

def parse_TCRpMHCI_pdbs():
    pass

def parse_TCRpMHCII_pdbs():
    pass
