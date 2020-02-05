#!/usr/bin/python
###    ###
from modelling_scripts import data_prep
import pickle
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
import glob
import time
from Bio.SubsMat import MatrixInfo
from random import choice

## TODO
## Select template structure by peptide BLAST positive score with PAM30
## make cmd_modeller output to one specific directory
## delete all useless produced files


### Retriving Dictionary with PDB IDs and chain lengths ###

##IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/final_mhc1_3d_structure_data_with_pdb_ids.tsv')
#IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/mhc_ligand_table_export_1578913026.csv', 76, 80, ',')

IDD_file = open('data/IDs_ChainsCounts_dict.pkl', 'rb')
IDD = pickle.load(IDD_file)
bad_IDs = pickle.load(IDD_file)
IDD_file.close()

### Organizing Allele IDs in a dictionary ###
allele_ID = {}
for key in IDD:
    if 'HLA' in IDD[key]['allele']:
        try:
            allele_ID[IDD[key]['allele']].append(key)
        except KeyError:
            allele_ID[IDD[key]['allele']] = [key]

pept_seqs = data_prep.get_peptides_from_csv('data/table_1_AnAnalysisofNaturalTCellResponsestoPredictedTumorNeoepitopes.csv', 1, 4, ',')

## for pept in pept_seqs: ### TODO later, parallelized on Cartesius

counts_dict = {'A1&A2' : 0, 'A1' : 0, 'A2' : 0, 'pos' : 0}
anch_dict = { 8: (1, 7), 9 : (1, 8), 10 : (1, 9), 
              11 : (1, 10), 12 : (1, 10)}
cutoff = 5

maxs = []
maxsl = []
print_results = False

for k, pept_seq in enumerate(pept_seqs):
    query = 'query_' + str(k + 1)
    
    pept = pept_seq[0]
    length = len(pept)
    #print('####################################################')
    #print('')
    #print('Peptide:  ', pept)
    #print('')
    
    
    scores_rules = { (anch_dict[length][0]-1) : {'pos' : 1, 'ID' : 2}, (anch_dict[length][0]+1) : {'pos' : 1, 'ID' : 2}, 
                    anch_dict[length][0] : {'pos' : 5, 'ID': 10}, anch_dict[length][1] : {'pos' : 3, 'ID' : 8}}
          
    max_pos = 0
    count = 0
    HLA_count = 0
    len_count = 0
    pos_list = []
    
    anch_1, anch_2 = anch_dict[length]
    for ID in IDD:
        count += 1
        score = 0
        temp_pept = IDD[ID]['pept_seq']
        if 'HLA' in IDD[ID]['allele']:                      ## Same Allele
            HLA_count += 1
            if length == len(temp_pept):                 ## Same Length
                len_count += 1
                for i, (aa, bb) in enumerate(zip(pept, temp_pept)):
                    ID_flag = False
                    pos_flag = False
                    try:
                        if aa == bb:
                            #score += 1 #
                            ID_flag = True
                        elif MatrixInfo.pam30[aa, bb] > 0:    #pos is identity plus positives
                            #score += 1 #
                            pos_flag = True
                        else:
                            pass
                    except KeyError:
                        try:
                            if MatrixInfo.pam30[bb, aa] > 0:
                                #score += 1 #
                                pos_flag = True
                            else:
                                pass
                        except KeyError:
                            #print('Broken peptide in IDD. Passing over.')
                            score = -50
                            pass
                    
                    if ID_flag == True:
                        try:
                            score += scores_rules[i]['ID']
                        except KeyError:
                            score += 1
                    elif pos_flag == True:
                        try:
                            score += scores_rules[i]['pos']
                        except KeyError:
                            score += 1
                    else:
                        pass

        if score > max_pos:
            max_pos = score
        if score > 0:
            pos_list.append((score, temp_pept, ID))
            
    max_list = []
    for pos in pos_list:
        if pos[0] == max_pos:
            max_list.append(pos)
            #print(pos)
    maxs.append(max_pos)
    maxsl.append(len(max_list))

    if len(max_list) == 1:
        template_ID = max_list[0][2]
    else:
        template_ID = choice(max_list)[2]
    if print_results:
        print('####################################################')
        print('')
        print('Peptide:  ', pept)
        print('')
        print('Pos:  ', pos)
        print('')
        print('####################################################')
    '''
    poss_flag = False
    A1_flag = False
    A2_flag = False
    A1A2_flag = False
    for pos in pos_list:
        if pos[0] == max_pos:
            max_list.append(pos)
            print(pos)
            if pos[1][anch_1] == pept[anch_1] and pos[1][anch_2] == pept[anch_2]:
                #print('Identical anchors ', pos)
                A1A2_flag = True
            elif pos[1][anch_1] == pept[anch_1] and pos[1][anch_2] != pept[anch_2]:
                #print('Identical P%i ' %anch_1, pos)
                A1_flag = True
            elif pos[1][anch_1] != pept[anch_1] and pos[1][anch_2] == pept[anch_2]:
                #print('Identical P%i ' %anch_2, pos)
                A2_flag = True
            elif pos[1][anch_1] != pept[anch_1] and pos[1][anch_2] == pept[anch_2]:
                poss_flag = True
                #print('Simply highest pos ', pos)
    if A1A2_flag == True:
        counts_dict['A1&A2'] += 1
    elif A1_flag == True:
        counts_dict['A1'] += 1
    elif A2_flag == True:
        counts_dict['A2'] += 1
    elif poss_flag == True:
        counts_dict['pos'] += 1
    '''
    
            
#%%
    ### Delete this Exception only when you are able to select the best template for each peptide
    #raise Exception("Pipeline in development. Don't run this block yet, it will only produce useless files")
    
    '''
    
    
    '''
    
    
    
    ### Producing a Fasta file for each putative template-target alignment ###
    sequences, empty_seqs = data_prep.get_pdb_seq([template_ID])
    
    template_seqr = SeqRecord(Seq(sequences[0]['M'], IUPAC.protein), id=template_ID, name = 'HLA' + ID)
    target_seqr = SeqRecord(Seq(sequences[0]['M'], IUPAC.protein), id=query, name = 'HLA_' + query)
    SeqIO.write((template_seqr, target_seqr), "data/FASTAs/%s.fasta" %template_ID, "fasta")
    peptide_seq = sequences[0]['P']
    
    final_alifile_name = 'data/Alignments/%s.ali' %template_ID
    
    final_alifile = open(final_alifile_name, 'w')
    i = 0
    for line in open('data/FASTAs/%s.fasta' %template_ID, 'r'):
        #print(line)
        if line.startswith('>') and i == 0:
            final_alifile.write('>P1;' + line.split(' ')[0].strip('>') + '\n')
            final_alifile.write('structure:data/PDBs/%s_MP_reres.pdb:1:M:9:P::::\n' %template_ID)
            i += 1
        elif line.startswith('>') and i == 1:
            final_alifile.write('/' + peptide_seq + '*')
            final_alifile.write('\n')
            final_alifile.write('\n>P1;' + line.split(':')[0].strip('>'))
            final_alifile.write('sequence:::::::::\n')
        else:
            final_alifile.write(line.rstrip())
    final_alifile.write('/' + str(pept) + '*')
    final_alifile.close()
    
    #os.system('rm data/FASTAs/%s.fasta' %template_ID)
    
    os.system('pdb_reres -1 data/PDBs/%s_MP.pdb > data/PDBs/%s_MP_reres.pdb' %(template_ID, template_ID))  # Renumbering the residues
    
    # Calculating all Atom contacts
    os.popen('modelling_scripts/contact-chainID_allAtoms data/PDBs/%s_MP_reres.pdb %s > data/all_contacts_%s.list' %(template_ID, cutoff, template_ID)).read()
    
    #Selecting only the anchors contacts
    anch_1_in = (anch_1 +1)
    anchor_1 = int(anch_1_in) + len(sequences[0]['M'])
    anch_2_in = (anch_2 +1)
    anchor_2 = anch_2_in +  len(sequences[0]['M'])
    
    real_anchor_2 = None
    anch_1_same = False
    anch_2_same = False
    
    if peptide_seq[anch_1_in-1] == pept[anch_1_in-1]:
        anch_1_same = True
    if anch_2_in == len(pept):
        if peptide_seq[anch_2_in-1] == pept[anch_2_in-1]:
            anch_2_same = True
    else:
        if peptide_seq[-1] == pept[anch_2_in-1]:
            anch_2_same = True
    
    if anchor_2 > (len(str(peptide_seq)) + len(peptide_seq)):
        real_anchor_2 = (len(str(peptide_seq)) + len(peptide_seq))
    
    ### Writing anchors contact list ###
    
    with open( 'data/all_contacts_%s.list' %template_ID, 'r') as contacts:
        with open('data/contacts_%s.list' %template_ID, 'w') as output:
            if real_anchor_2:                                                                     ### If the target peptide is longer than the templtate peptide
                for line in contacts:
                    #print(line[30:33])
                    p_aa_id = line.split("\t")[7]                                                 ### position id of the template peptide residue
                    p_atom = line.split("\t")[8]                                                  ### atom name of the template peptide residue
                    m_aa_id = (line.split("\t")[2]).split(' ')[0]
                    if anch_1_same == True:                                                       ### If the target anchor 1 residue is the same as the template anchor 1 residue
                        if int(p_aa_id) == anchor_1:
                            output.write(line)
                    else:
                        if int(p_aa_id) == anchor_1 and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line)
                    if anch_2_same == True:                                                       ### If the target anchor 2 residue is the same as the template anchor 2 residue
                        if int(p_aa_id) == real_anchor_2:
                            output.write(line[:30] + str(anchor_2) + line[34:])
                    else:
                        if int(p_aa_id) == real_anchor_2 and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line[:30] + str(anchor_2) + line[34:])
            else:
                for line in contacts:
                    #print(line.split("\t"))
                    p_aa_id = line.split("\t")[7]
                    p_atom = line.split("\t")[8]
                    m_aa_id = (line.split("\t")[2]).split(' ')[0]
                    if anch_1_same == True:                                                       ### If the target anchor 1 residue is the same as the template anchor 1 residue
                        if int(p_aa_id) == anchor_1:
                            output.write(line)
                    else:
                        if int(p_aa_id) == anchor_1 and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line)
                    if anch_2_same == True:                                                       ### If the target anchor 2 residue is the same as the template anchor 2 residue
                        if int(p_aa_id) == anchor_2:
                            output.write(line)
                    else:
                        if int(p_aa_id) == anchor_2 and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line)
                    #if (int(p_aa_id) == anchor_1 or int(p_aa_id) == anchor_2) and ('CA' in p_atom or 'CB' in p_atom):
                    #        output.write(line)
    
    
    with open('data/instructions.txt', 'w') as instr_file:
        instr_file.write(template_ID + '\n')
        instr_file.write(str(anchor_1) + '\n')
        instr_file.write(str(anchor_2))
    
    #Finally launching Modeller. Hopefully.
    
    #command = ['python', 'modelling_scripts/cmd_modeller.py', final_alifile_name, template_ID, '>3ROO:A']
    
    #proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
    #(out, err) = proc.communicate()
    
    os.popen('/usr/bin/python2.7 modelling_scripts/cmd_modeller.py %s %s %s' %(final_alifile_name, template_ID, query)).read()

    raise Exception('OK.')
#os.popen('/usr/bin/python2.7 modelling_scripts/cmd_modeller.py %s 1k5n_MP query_1ogt_MP' %final_alifile_name).read()
