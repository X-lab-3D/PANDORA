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
from Bio.PDB import PDBParser
from Bio.Data.SCOPData import protein_letters_3to1 as to_one_letter_code
from modelling_scripts import del_uncommon_residues_pdbs as durp

### OPEN CSV FILE WITH PDB CODES ###
#def pdb_reres_allchains(pdb_file)


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

def download_ids_alleles_imgt(out_tsv = 'auto_generated_IDs_alleles_from_IMGT.tsv', out_pkl = 'IDs_and_alleles_identity_percs_from_imgt.pkl', print_outfiles = False):
    cwd = os.getcwd()
    
    '''
    params = { 'ReceptorType' : 'peptide/MH1',
            'type-entry': 'PDB'}
    '''
    params = { 'ReceptorType' : 'MH1',
        'type-entry': 'PDB'}
    
    url = "http://www.imgt.org/3Dstructure-DB/cgi/3Dquery.cgi"
    
    data = urllib.parse.urlencode(params)
    data = data.encode('ascii') # data should be bytes
    req = urllib.request.Request(url, data)
    
    IDs_list = []
    
    
    with urllib.request.urlopen(req) as response:
        text = response.read().decode('utf-8')
        text = text.splitlines()
        if print_outfiles: 
            temp_outfile = open('../outputs/test/test_download/test.html', 'w')
            temp_outfile.write(text)
            temp_outfile.close()
    
    IDs_list = [x for x in text if 'href' in x and 'pdbcode' in x]
    IDs_list = [x.split('"') for x in IDs_list]
    IDs_list = [x[3][-4:] for x in IDs_list]
    

    ID_allele = {}
    for ID in IDs_list:
        
        str_url = 'http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode=%s' %ID
        #str_url = "http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode=6ujo"
        
        req = urllib.request.Request(str_url)
        
        with urllib.request.urlopen(req) as response: #TODO: Handle exception "HTTPError: Internal Server Error"
            text = response.read().decode('utf-8')
            text = text.splitlines()
            if print_outfiles: 
                temp_outfile = open('%s/outputs/test/test_download/%s.html' %(cwd,ID), 'w')
                temp_outfile.write(text)
                temp_outfile.close()
        
        domain_flag = False
        domains = []
        domain = []
        for line in text:
            if 'collier_perles' in line:
                continue
            if line == '<td class="titre_h" align="center" rowspan="6">G-DOMAIN</td>':
                domain_flag = True
                #print(line)
            elif line =='\n' or line =='':
                domain_flag = False
                if domain != []:
                    domains.append(domain)
                    domain = []
            if domain_flag:
                domain.append(line)
        
        for i, dom in enumerate(domains):
            flag = False
            descr = 0
            for j, data in enumerate(dom):
                if flag == True and j == (descr + 2):
                    #print('Is this what are you looking for? :', data)
                    flag = False
                    domains[i] = data.split('title')[0]
                    continue
                if 'IMGT gene and allele name' in data:
                    descr = j
                    flag = True
        
        percs = [{}, {}]
        clean_domains = deepcopy(domains)
        for i, domain in enumerate(domains):
            clean_domains[i] = domains[i].replace(',','').replace(';','').split('&nbsp')
        for i, domain in enumerate(clean_domains[:2]):
            for j, d in enumerate(clean_domains[i]):
                if '%' in d:
                    percs[i][clean_domains[i][j-1]] = float((nestedExpr('(',')').parseString(d).asList())[0][0].replace('%',''))
            clean_domains[i] = [x for x in clean_domains[i] if '%' not in x]
        
        
        # TODO: in future: in case of ambiguity, assign allele to the templates according to patient genotype?
        
        if percs != [{}, {}]:
            set_percs = set(percs[0]).intersection(*percs)
            
                
        ID_allele[ID] = [set_percs, percs] #TODO: check if allele is in both domains
        
    if out_pkl:
        with open('%s/data/csv_pkl_files/%s' %(cwd, out_pkl), 'wb') as outpkl:
            pickle.dump(ID_allele, outpkl)
    
    if out_tsv:
        with open('%s/data/csv_pkl_files/%s' %(cwd, out_tsv), 'w') as outtsv:
            outtsv.write('PDB ID' + '\t' + 'Equal identity alleles' + '\n')
            for ID in ID_allele:
                if len(ID_allele[ID][0]) != 0:
                    outtsv.write(ID + '\t' + (';').join([x for x in ID_allele[ID][0]]) +'\n')
                elif len(ID_allele[ID][0]) == 0:
                    to_write = []
                    for x in ID_allele[ID][1]:
                        to_write += list(x.keys())
                    to_write = list(set(to_write))
                    outtsv.write(ID + '\t' + (';').join(to_write) +'\n')

def imgt_retrieve_clean(ids_filename, id_clmn, allele_clmn, delimiter, empty_rows = [0]):


    IDs = []
    with open(ids_filename, 'r') as idsfile:
        spamreader = csv.reader(idsfile, delimiter=delimiter)
        for i, row in enumerate(spamreader):
            if i in empty_rows:
                pass
            else:
                ID = row[id_clmn]
                allele = row[allele_clmn].split(';')
                IDs.append((ID, allele))
    cwd = os.getcwd()

    if "contact-chainID_allAtoms" not in os.listdir('%s/modelling_scripts' %cwd):
        os.popen('g++ %s/modelling_scripts/contact-chainID_allAtoms.cpp -o %s/modelling_scripts/contact-chainID_allAtoms' %(cwd, cwd)).read()
    #print('CWD:', cwd)
    bad_IDs = {}
    err_1 = 0
    err_2 = 0
    err_3 = 0
    err_4 = 0
    err_5 = 0
    err_6 = 0
    err_7 = 0
    IDs_dict = {}
    print('Started fetching URLs')
    for entry in IDs:
        ID = entry[0]
        ID = ID.upper()
        ID = ID.rstrip()
        #if '%s_MP.pdb' %ID in os.listdir('%s/data/PDBs/' %cwd):
        #    continue
        url = 'http://www.imgt.org/3Dstructure-DB/IMGT-FILE/IMGT-%s.pdb.gz' %ID
        filepath = '%s/data/PDBs/%s.pdb' %(cwd, ID)
        print('Fetching ', url)

        try:
            urllib.request.urlretrieve( url, filepath + '.gz')
        except: #urllib.error.HTTPError:
            bad_IDs[ID] = '#1'
            print('ERROR TYPE #1: URL not found. %s added to Bad_IDs' %ID)
            print('##################################')
            err_1 += 1
            continue
        with gzip.open('data/PDBs/%s.pdb.gz' %ID, 'rb') as f_in:
            with open('data/PDBs/%s.pdb' %ID, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        subprocess.check_call(['rm', 'data/PDBs/%s.pdb.gz' %ID])
        selected_chains_filepath = '%s/data/PDBs/%s_AC.pdb'%(cwd, ID)
        a_renamed_filepath = '%s/data/PDBs/%s_MC.pdb'%(cwd, ID)
        new_filepath = '%s/data/PDBs/%s_MP.pdb'%(cwd, ID)


        #print('Opening %s' %ID)
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
        except (IndexError, ValueError):
            bad_IDs[ID] = '#2'
            print('ERROR TYPE #2: Something went wrong with PDBParser. %s added to Bad_IDs' %ID)
            print('##################################')
            err_2 += 1
            subprocess.check_call(['rm', 'data/PDBs/%s.pdb' %ID])
            continue

        ### Selecting Alpha chain and peptide chain

        putative_pID = {}
        for chain in seqs:
            if type(seqs[chain]) == str:
                length = len(seqs[chain])
                if chain != ' ':               #This might be removed, since get_seqs() is  not outputting anymore with empty chain ids
                    count_dict[chain] = length
                if length > 250 and length < 300 and not aID:
                    aID = chain
                    #print('aID : ', aID)
                    #print('length : ', length)
                elif length > 7 and length < 25:
                    putative_pID[chain] = 0
                    #print('pID : ', pID)
                    #print('length : ', length)
        if len(putative_pID) == 1:
            pID = list(putative_pID.keys())[0]
        elif len(putative_pID) > 1:
            os.system('./modelling_scripts/contact-chainID_allAtoms %s 5 > ./data/dist_files/%s.dist' %(pdbf, ID))
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


        try:
            if aID.islower() or pID.islower():
                bad_IDs[ID] = '#3'
                print('ERROR TYPE #3: Lower case chain ID. %s added to Bad_IDs' %ID)
                print('##################################')
                err_3 += 1
                subprocess.check_call(['rm', 'data/PDBs/%s.pdb' %ID])
                continue
            else:
                pass
        except AttributeError:
            bad_IDs[ID] = '#4'
            print('ERROR TYPE #4: False in chain ID, one chain does not respect settled lenght criteria. %s added to Bad_IDs' %ID)
            print(count_dict)
            print('##################################')
            err_4 += 1
            subprocess.check_call(['rm', 'data/PDBs/%s.pdb' %ID])
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
                bad_IDs[ID] = '#6'
                print('ERROR TYPE #6: Unrecognized chain lenght. %s added to Bad_IDs' %ID)
                print('##################################')
                err_6 += 1
                subprocess.check_call(['rm', 'data/PDBs/%s.pdb' %ID])
                continue
        if Atcr != 0 or Btcr != 0:
            bad_IDs[ID] = '#5'
            print('ERROR TYPE #5: TCRs in structure. %s added to Bad_IDs' %ID)
            print('##################################')
            err_5 += 1
            subprocess.check_call(['rm', 'data/PDBs/%s.pdb' %ID])
            continue
        ######
        else:
            count_dict['allele'] = entry[1]
            count_dict['pept_seq'] = seqs[pID]
            IDs_dict[ID] = count_dict
            os.system('pdb_selchain -%s,%s %s > %s' %(aID, pID, filepath, selected_chains_filepath))
            if pID == 'M':
                os.system('pdb_rplchain -%s:P %s > %s' %(pID, selected_chains_filepath, a_renamed_filepath))
                os.system('pdb_rplchain -%s:M %s > %s' %(aID, a_renamed_filepath, new_filepath))
            else:
                os.system('pdb_rplchain -%s:M %s > %s' %(aID, selected_chains_filepath, a_renamed_filepath))
                os.system('pdb_rplchain -%s:P %s > %s' %(pID, a_renamed_filepath, new_filepath))
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

    ### ADD CORRECTION FOR SEP, F2F, etc.
    os.system('python ./tools/change_sep_in_ser.py ./data/PDBs/')
    print('Removing uncommon residue files')
    uncommon_pdbs = durp.move_uncommon_pdbf('./data/PDBs/')
    for u_pdb in uncommon_pdbs:
        try:
            del IDs_dict[u_pdb]
            bad_IDs[u_pdb] = '#7'
            err_7 += 1
        except KeyError:
            print("Tried to delete %s from dataset error #7, error encountered" %u_pdb)

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
    remove_error_templates()
    return IDs_dict, bad_IDs

def remove_error_templates():
    try:
        os.system('mkdir data/PDBs//unused_templates')
    except:
        pass
    with open('data/csv_pkl_files/IDs_ChainsCounts_dict.pkl', 'rb') as inpkl:
        IDD = pickle.load(inpkl)
        bad_IDs = pickle.load(inpkl)
        inpkl.close()
    with open('data/csv_pkl_files/error_templates.tsv') as infile: #To be removed later
        r = csv.reader(infile, delimiter='\t')
        for i, row in enumerate(r):
            if i != 0:
                try:
                    del IDD[row[0]]
                    bad_IDs[row[0]] = row[1]
                    os.system('mv data/PDBs/pMHCI/%s_MP.pdb data/unused_templates/' %row[0])
                except:
                    pass
        infile.close()
    with open("data/csv_pkl_files/IDs_ChainsCounts_dict.pkl", "wb") as outpkl:
        pickle.dump(IDD, outpkl)
        pickle.dump(bad_IDs, outpkl)
        outpkl.close()
        return IDD, bad_IDs


def get_pdb_seq(IDs):
    sequences = []
    empty_seqs = []
    for ID in IDs:
        pdbf = 'data/PDBs/pMHCI/%s_MP.pdb' %ID
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
