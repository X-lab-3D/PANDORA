#!/usr/bin/python
###    ###

import os
import sys
import time
import copy
import pickle
import csv
#import subprocess

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SubsMat import MatrixInfo
#from Bio.Align import substitution_matrices
#PAM30 = substitution_matrices.load('PAM30')


from random import choice
from joblib import Parallel, delayed
from multiprocessing import Manager

#PANDORA modules
sys.path.append('/home/dariom/PANDORA_master_to_package/') #TODO: change in final release
import PANDORA
from PANDORA.parsing import utils

def check_model_existance(k, outdir_name, filename_start, filename_end):
    '''
    Checks if models have been already generated for the current run and case
    '''
    ext_flag = False
    for folder in os.listdir('PANDORA_files/outputs/%s' %outdir_name):
        if 'query' in folder:
            if str(k+1) == folder.split('_')[2]:
                print('WARNING: Existing directory for this query.')
                for name in os.listdir('PANDORA_files/outputs/%s/%s' %(outdir_name, folder)):
                    if filename_start in name and filename_end in name:
                        print('WARNING: This query has already been modelled. Moving on.')
                        ext_flag = True
                        break
                if ext_flag == True:
                    break

    if ext_flag == True:
        print('Exiting here.')
        return True
    else:
        return False

def allele_name_adapter(allele, allele_ID):
    '''
    Cuts the given allele name to make it consistent with the alleles in allele_ID.

    Args:
        allele(str) : Allele name
        allele_ID(dict) : Dictionary of structure IDs (values) in the dataset for each allele (keys)
    '''
    #homolog_allele = '--NONE--'
    for a in range(len(allele)):
        if allele[a].startswith('HLA'):      # Human
            if any(allele[a] in key for key in list(allele_ID.keys())):
                pass
            elif any(allele[a][:8] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:8]
            elif any(allele[a][:6] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:6]
            else:
                allele[a] = allele[a][:4]
        elif allele[a].startswith('H2'):    # Mouse
            #homolog_allele = 'RT1'
            if any(allele[a] in key for key in list(allele_ID.keys())):
                pass
            elif any(allele[a][:4] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:4]
            else:
                allele[a] = allele[a][:3]
        elif allele[a].startswith('RT1'):          # Rat
            #homolog_allele = 'H2'
            if any(allele[a] in key for key in list(allele_ID.keys())):
                pass
            elif any(allele[a][:5] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:5]
            else:
                allele[a] = allele[a][:4]
        elif allele[a].startswith('BoLA'):        # Bovine
            if any(allele[a] in key for key in list(allele_ID.keys())):
                pass
            elif any(allele[a][:10] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:10]
            elif any(allele[a][:7] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:7]
            else:
                allele[a] = allele[a][:5]
        elif allele[a].startswith('SLA'):        # Suine
            if any(allele[a] in key for key in list(allele_ID.keys())):
                pass
            elif any(allele[a][:9] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:9]
            elif any(allele[a][:6] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:6]
            else:
                allele[a] = allele[a][:4]
        elif allele[a].startswith('MH1-B'):        # Chicken
            if any(allele[a] in key for key in list(allele_ID.keys())):
                pass
            elif any(allele[a][:8] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:8]
            else:
                allele[a] = allele[a][:6]
        elif allele[a].startswith('BF2'):        # Chicken
            if any(allele[a] in key for key in list(allele_ID.keys())):
                pass
            elif any(allele[a][:6] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:6]
            else:
                allele[a] = allele[a][:4]
        elif allele[a].startswith('Mamu'):       # Monkey
            if any(allele[a] in key for key in list(allele_ID.keys())):
                pass
            elif any(allele[a][:13] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:13]
            elif any(allele[a][:9] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:9]
            else:
                allele[a] = allele[a][:5]
        elif allele[a].startswith('Eqca'):        # Horse
            if any(allele[a] in key for key in list(allele_ID.keys())):
                pass
            elif any(allele[a][:10] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:10]
            elif any(allele[a][:7] in key for key in list(allele_ID.keys())):
                allele[a] = allele[a][:7]
            else:
                allele[a] = allele[a][:5]
    return(allele)#, homolog_allele)

def select_template(IDD, target_id, allele, length, pept, print_results): #homolog_allele,
    '''
    Select template for p:MHC I Modelling.

    Args:
        IDD(dict) : dictionary containing all the templates ID and relative informations
        allele(str) : target allele name
        length(int) : target peptide length
        pept(str) : target peptide sequence
        print_result(bool) : if True prints out the selected template structure and peptide

    '''

    putative_templates = []
    max_pos = -1000
    pos_list = []

    for ID in IDD:
        for a in allele:
            if any(a in key for key in IDD[ID]['allele']): #or homolog_allele in IDD[ID]['allele']:                       ## Same Allele
                putative_templates.append(ID)
    putative_templates = list(set(putative_templates))


    for ID in putative_templates:
        score = 0
        temp_pept = IDD[ID]['pept_seq']
        min_len = min([length, len(temp_pept)])
        score -= ((abs(length - len(temp_pept)) ** 2.4)) #!!!  ## Gap Penalty
        for i, (aa, bb) in enumerate(zip(pept[:min_len], temp_pept[:min_len])):
            try:
                gain = MatrixInfo.pam30[aa, bb]
                score += gain
            except KeyError:
                try:
                    gain = MatrixInfo.pam30[bb, aa]
                    score += gain
                except KeyError:
                    score = -50
                    pass

        if score > max_pos:
            max_pos = score
        pos_list.append((score, temp_pept, ID))

    max_list = []
    for pos in pos_list:
        if pos[0] == max_pos:
            max_list.append(pos)
    #maxs.append(max_pos)
    #maxsl.append(len(max_list))

    if len(max_list) == 0:
        return (target_id, pept, allele, "NA", "NA", "NA", "No positive scoring template peptides")
    elif len(max_list) == 1:
        template = max_list[0]
        template_ID = template[2]
        template_pept = template[1]
    else:
        template = (choice(max_list))
        template_ID = template[2]
        template_pept = template[1]
    if print_results:
        print('####################################################')
        print('')
        print('Peptide:  ', template[1])
        print('')
        print('Pos:  ', template)
        print('')
        print('####################################################')

    return((template, template_ID, template_pept))

def make_ali_files(templ_sequences, target_sequence, template_ID, target_id, 
                   query, template_pept, temp_anch_1, temp_anch_2, pept, 
                   anch_1, anch_2, length, remove_temp_outputs=True):
    '''
    Aligns target sequence and template sequence
    '''

    template_seqr = SeqRecord(Seq(templ_sequences['M']), id=template_ID, name = template_ID)
    target_seqr = SeqRecord(Seq(target_sequence), id=query, name = query)
    SeqIO.write((template_seqr, target_seqr), "%s.fasta" %template_ID, "fasta")

    os.system('muscle -in %s.fasta -out %s.afa -quiet' %(template_ID, template_ID))
    os.system('rm %s.fasta' %template_ID)

    ali_template_pept, ali_target_pept = utils.align_peptides(template_pept, temp_anch_1, temp_anch_2, pept, anch_1, anch_2)

    final_alifile_name = '%s.ali' %template_ID
    final_alifile = open(final_alifile_name, 'w')
    i = 0
    target_ini_flag = False
    template_ini_flag = False
    modeller_renum = 1
    for line in open('%s.afa' %template_ID, 'r'):
        if line.startswith('>') and i == 0:
            final_alifile.write('>P1;' + line.split(' ')[0].strip('>') + '\n')
            final_alifile.write('structure:%s/PDBs/pMHCI/%s_MP.pdb:%s:M:%s:P::::\n' %(PANDORA.PANDORA_data, template_ID, str(templ_sequences['M_st_ID']), str(len(ali_template_pept))))
            template_ini_flag = True
            i += 1
        elif line.startswith('>') and i == 1:
            final_alifile.write('/' + ali_template_pept + '*')
            final_alifile.write('\n')
            final_alifile.write('\n>P1;' + line.split(':')[0].strip('>'))
            final_alifile.write('sequence:::::::::\n')
            target_ini_flag = True
        else:
            final_alifile.write(line.rstrip())
    final_alifile.write('/' + str(ali_target_pept) + '*')
    final_alifile.close()

    if remove_temp_outputs:
        os.system('rm *.afa')
    return(final_alifile_name)

def write_ini_scripts(anch_1, anch_2, template_ID, final_alifile_name, query):
        with open('MyLoop.py', 'w') as myloopscript:
            MyL_temp = open(PANDORA.PANDORA_path + '/modelling/MyLoop_template.py', 'r')
            for line in MyL_temp:
                if 'self.residue_range' in line:
                    myloopscript.write(line %(anch_1 + 2, anch_2))
                elif 'SPECIAL_RESTRAINTS_BREAK' in line:
                    break
                elif 'contact_file = open' in line:
                    myloopscript.write(line %template_ID)
                else:
                    myloopscript.write(line)
            MyL_temp.close()

        with open('cmd_modeller_ini.py', 'w') as modscript:
            cmd_m_temp = open(PANDORA.PANDORA_path + '/modelling/cmd_modeller_ini.py', 'r')
            for line in cmd_m_temp:
                if 'alnfile' in line:
                    modscript.write(line %final_alifile_name)
                elif 'knowns' in line:
                    modscript.write(line %(template_ID, query))
                else:
                    modscript.write(line)
            cmd_m_temp.close()

def write_contact_list(pept, anch_1, anch_2, template_pept, temp_anch_1, temp_anch_2, template_ID):
    
    ### Select only the anchors contacts
    anch_1_same = False
    anch_2_same = False

    if template_pept[temp_anch_1] == pept[anch_1]:
        anch_1_same = True
    if template_pept[temp_anch_2] == pept[anch_2]:
        anch_2_same = True
    if pept[anch_1] == 'G':
        first_gly = True
    else:
        first_gly= False
    if pept[anch_2] == 'G':
        last_gly = True
    else:
        last_gly = False

    with open( 'all_contacts_%s.list' %template_ID, 'r') as contacts:                             # Template contacts
        with open('contacts_%s.list' %template_ID, 'w') as output:
            for line in contacts:
                p_aa_id = line.split("\t")[7]                                     ### position id of the template peptide residue
                p_atom = line.split("\t")[8]                                      ### atom name of the template peptide residue
                if anch_1_same == True:                                           ### If the target anchor 1 residue is the same as the template anchor 1 residue
                    if int(p_aa_id) == (anch_1+1):
                        output.write(line)
                else:
                    if int(p_aa_id) == (anch_1+1):
                        if first_gly:
                            if 'CA' in p_atom:
                                output.write(line)
                        else:
                            if 'CA' in p_atom or 'CB' in p_atom:
                                output.write(line)
                if anch_2_same == True:                                                       ### If the target anchor 2 residue is the same as the template anchor 2 residue
                    if int(p_aa_id) == (anch_2+1):
                        output.write(line)
                else:
                    if int(p_aa_id) == (anch_2+1):
                        if last_gly:
                            if 'CA' in p_atom:
                                output.write(line)
                        else:
                            if 'CA' in p_atom or 'CB' in p_atom:
                                output.write(line)

def write_modelling_scripts(anch_1, anch_2, template_ID, target_id, final_alifile_name):
    with open('MyLoop.py', 'w') as myloopscript:
        MyL_temp = open(PANDORA.PANDORA_path + '/modelling/MyLoop_template.py', 'r')
        for line in MyL_temp:
            if 'self.residue_range' in line:
                myloopscript.write(line %(anch_1 + 2, anch_2))
            elif 'contact_file = open' in line:
                myloopscript.write(line %template_ID)
            else:
                myloopscript.write(line)
        MyL_temp.close()

    with open('cmd_modeller.py', 'w') as modscript:
        cmd_m_temp = open(PANDORA.PANDORA_path + '/modelling/cmd_modeller_template.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line %final_alifile_name)
            elif 'knowns' in line:
                modscript.write(line %(template_ID, target_id))
            else:
                modscript.write(line)
        cmd_m_temp.close()
