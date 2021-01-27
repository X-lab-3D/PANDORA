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
#from random import choice

###IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/csv_pkl_files/final_mhc1_3d_structure_data_with_pdb_ids.tsv')

#IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/csv_pkl_files/mhcI_structures_IEDB.csv', 42, 43, ';', [0,1])

### Retriving Dictionary with PDB IDs and chain lengths ###

IDD_file = open('data/csv_pkl_files/IDs_ChainsCounts_dict.pkl', 'rb')
IDD = pickle.load(IDD_file)
bad_IDs = pickle.load(IDD_file)
IDD_file.close()

remove_temp_outputs = False
cutoff = 5
outdir_name = 'table_1_v2'

### Organizing Allele IDs in a dictionary ###
allele_ID = {}
for key in IDD:
    try:
        allele_ID[IDD[key]['allele']].append(key)
    except KeyError:
        allele_ID[IDD[key]['allele']] = [key]

run_dict = {0: ('HLA-A*02:01', '1duz', '3MRB', '2', '9'), 1: ('HLA-A*02:01', '1tvb', '3MRB', '2', '9'), 
            2: ('HLA-A*02:01', '2v2w', '3MRB', '2', '9'), 3: ('HLA-A*02:01', '3pwn', '3MRB', '2', '9'),
            4: ('HLA-A*02:01', '4k7f', '3MRB', '2', '9'), 5: ('HLA-B*27:09', '1ogt', '1K5N', '2', '9'),
            6: ('HLA-B*27:09', '1w0w', '1K5N', '2', '9'), 7: ('HLA-B*27:09', '2a83', '1K5N', '2', '9'),
            8: ('HLA-B*27:09', '2bst', '1K5N', '2', '9'), 9: ('HLA-B*27:09', '3bp4', '1K5N', '2', '9'),
            10:('HLA-B*35:01', '1a9e', '2CIK', '2', '9'), 11:('HLA-B*35:01', '3lkn', '2CIK', '2', '9'),
            12:('HLA-B*35:01', '3lkp', '2CIK', '2', '9'), 13:('HLA-B*35:01', '3lkr', '2CIK', '2', '9'),
            14:('HLA-B*35:01', '3lks', '2CIK', '2', '9'), 15:('HLA-B*44:02', '1m6o', '3KPM', '2', '9'),
            16:('HLA-B*44:02', '1sys', '3KPM', '2', '9'), 17:('HLA-B*44:02', '3kpp', '3KPM', '2', '9'),
            18:('HLA-B*44:02', '3l3j', '3KPM', '2', '9'), 19:('HLA-B*44:02', '3l3k', '3KPM', '2', '9'),
            20:('HLA-E*01:01', '1kpr', '3BZF', '2', '9'), 21:('HLA-E*01:01', '1ktl', '3BZF', '2', '9'),
            22:('HLA-E*01:01', '2esv', '3BZF', '2', '9'), 23:('HLA-E*01:01', '3am8', '3BZF', '2', '9'),
            24:('HLA-E*01:01', '3bze', '3BZF', '2', '9')}
for run in range(25):
    if run != 0:
        os.chdir('../../../')
    
    ### Retriving target allele and sequence file ###
    '''
    print('##############################')
    print('Please enter your allele name. You can choose between the ones listed below')
    print('##############################')
    print('')
    print('')
    print(allele_ID.keys())
    '''
    #allele = input()
    allele = run_dict[run][0]
    #allele, inseq_file = 'FLA-E*01801', 'data/5xmf.fasta'
    #allele, inseq_file = 'HLA-B*27:09', 'data/1ogt.fasta'
    
    '''
    print('##############################')
    print('Please select your target sequence file.')
    print('##############################')
    print('')
    
    for file in glob.glob(os.getcwd() + '/data/Targets/*.fasta'):
        print(file.split('/')[-1].split('.')[0])
    '''
    #inseq_file = 'data/5xmf.fasta'
    #inseq_file = 'data/1tom.fasta'
    #target_name = input()
    target_name = run_dict[run][1]
    inseq_file = 'data/Targets/%s.fasta' %target_name
    inseqs = list(SeqIO.parse(inseq_file, 'fasta'))
    A_target = inseqs[0]
    P_target = inseqs[1]
    
    ### Producing a Fasta file for each putative template-target alignment ###
    sequences, empty_seqs = data_prep.get_pdb_seq(allele_ID[allele])
    for bad_ID in empty_seqs:
        del IDD[bad_ID]
        allele_ID[allele].remove(bad_ID)
    
    ID_seqs_dict = {}
    
    ### Identifying the template ID and its peptide sequence ###
    if len(allele_ID[allele]) == 1:   ### In case we have only one template for this allele ###
        template_ID = allele_ID[allele][0]
        outdir = ('outputs/%s/%s_%s' %(outdir_name, template_ID.lower(), target_name))
        try:
            os.mkdir(outdir)
        except FileExistsError:
            print('WARNING: You are writing into an existing directory.')
            pass
        os.chdir(outdir)
    
        for i, ID in enumerate(allele_ID[allele]):
            #print(i, ID, sequences[i])
            ID_seqs_dict[ID] = ((sequences[i]['P'], len(sequences[i]['M'])))
            if ID == template_ID:
                template = SeqRecord(Seq(sequences[i]['M'], IUPAC.protein), id=ID, name = allele + ID)
                SeqIO.write((template, A_target), "%s.fasta" %ID, "fasta")
        template_pept = ID_seqs_dict[template_ID][0]
        os.system('muscle -in %s.fasta -out %s.afa' %(template_ID, template_ID))
        os.system('rm %s.fasta' %template_ID)
    
    else:                             ### In case we have multiple templates for this allele ###
        score_dict = {}
        max_id = 0
        templates_dict = allele_ID[allele]
        '''
        print('##############################')
        print('Please select one template between the ones listed below. You can enter the key number or the ID, if you are sure it is present')
        print('##############################')
        print('')
        print('')
        for i, j in enumerate(allele_ID[allele]):
            print(i,j)
        '''
        #temp_in = input()
        temp_in = run_dict[run][2]
        try:
            template_ID = allele_ID[allele][int(temp_in)]
        except:
            template_ID = temp_in
        outdir = ('outputs/%s/%s_%s' %(outdir_name, template_ID.lower(), target_name))
        try:
            os.mkdir(outdir)
        except FileExistsError:
            print('WARNING: You are writing in an existing directory.')
            pass
        os.chdir(outdir)
    
        for i, ID in enumerate(allele_ID[allele]):
            #print(i, ID, sequences[i])
            ID_seqs_dict[ID] = ((sequences[i]['P'], len(sequences[i]['M'])))
            if ID == template_ID:
                template = SeqRecord(Seq(sequences[i]['M'], IUPAC.protein), id=ID, name = allele + ID)
                SeqIO.write((template, A_target), "%s.fasta" %ID, "fasta")
    
        alifile = '%s.afa' %template_ID
        muscle_commands = ['muscle', '-in', '%s.fasta' %template_ID, '-out', '%s' %alifile, '-quiet']
        subprocess.check_call(muscle_commands)
        
        if remove_temp_outputs:
            os.system('rm %s.fasta' %template_ID)
    
        template_pept = ID_seqs_dict[template_ID][0]
    
    ### Writing a final .ali file with Template sequence / template pept sequence ; Target sequence / target pept sequence ###
    
    final_alifile_name = '%s.ali' %template_ID
    final_alifile = open(final_alifile_name, 'w')
    i = 0
    target_ini_flag = False
    template_ini_flag = False
    
    for line in open('%s.afa' %template_ID, 'r'):
        #print(line)
        if line.startswith('>') and i == 0:
            final_alifile.write('>P1;' + line.split(' ')[0].strip('>') + '\n')
            final_alifile.write('structure:../../../data/PDBs/%s_MP.pdb:%s:M:%s:P::::\n' %(template_ID, str(sequences[0]['M_st_ID']), str(len(template_pept)))) 
            #('structure:../../data/PDBs/%s_MP.pdb:1:M:9:P::::\n' %template_ID)
            template_ini_flag = True
            i += 1
        elif line.startswith('>') and i == 1:
            final_alifile.write('/' + template_pept + '*')
            final_alifile.write('\n')
            final_alifile.write('\n>P1;' + line.split(':')[0].strip('>') + '\n')
            final_alifile.write('sequence:::::::::\n')
            target_ini_flag = True
        else:
            final_alifile.write(line.rstrip())
    
    final_alifile.write('/' + str(P_target.seq) + '*')
    final_alifile.close()
    
    if remove_temp_outputs:
        os.system('rm *.afa')
    '''
    with open('instructions.txt', 'w') as instr_file:
        instr_file.write(template_ID + ' ' + str(1) + ' ' + str(9))
    os.popen('python2.6 ../../modelling_scripts/cmd_modeller_ini.py %s %s %s' %(final_alifile_name, template_ID, target_name.upper())).read()
    '''
    
    #Selecting only the anchors contacts
    '''
    print('##############################')
    print('Please input, one per time, the anchor positions (only the int number)')
    print('##############################')
    print('')
    print('The peptide of the selected template is %i residue long. The target one is %i residue long' %(len(str(template_pept)), len(str(P_target.seq))))
    print('')
    print('The template peptide sequence is %s, the target peptide seuqnce is %s' %(str(template_pept), str(P_target.seq)))
    '''
    print('ANCHOR 1:')
    #anchor_1 = int(input())
    anchor_1 = int(run_dict[run][3])
    #anchor_1 = int(anch_1_in) + ID_seqs_dict[template_ID][1]
    print('ANCHOR 2:')
    #anchor_2 = int(input())
    anchor_2 = int(run_dict[run][4])
    #anchor_2 = anch_2_in + ID_seqs_dict[template_ID][1]
    
    real_anchor_2 = None
    anch_1_same = False
    anch_2_same = False
    
    if template_pept[anchor_1-1] == P_target.seq[anchor_1-1]:
        anch_1_same = True
    if anchor_2 == len(P_target.seq):
        if template_pept[anchor_2-1] == P_target.seq[anchor_2-1]:
            anch_2_same = True
    else:
        if template_pept[-1] == P_target.seq[anchor_2-1]:
            anch_2_same = True
    
    if anchor_2 > (len(str(template_pept)) + ID_seqs_dict[template_ID][1]):
        real_anchor_2 = (len(str(template_pept)) + ID_seqs_dict[template_ID][1])
    
    with open('MyLoop.py', 'w') as myloopscript:
        MyL_temp = open('../../../modelling_scripts/MyLoop_template.py', 'r')
        for line in MyL_temp:
            
            if 'self.residue_range' in line:
                myloopscript.write(line %(anchor_1 + 2, anchor_2))
            elif 'SPECIAL_RESTRAINTS_BREAK' in line:
                break
            elif 'contact_file = open' in line:
                myloopscript.write(line %template_ID)
            else:
                myloopscript.write(line)
        MyL_temp.close()
    
    with open('cmd_modeller_ini.py', 'w') as modscript:
        cmd_m_temp = open('../../../modelling_scripts/cmd_modeller_ini.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line %final_alifile_name)
            elif 'knowns' in line:
                modscript.write(line %(template_ID, target_name.upper()))
            else:
                modscript.write(line)
        cmd_m_temp.close()
    
    os.popen('python3 cmd_modeller_ini.py').read()
    
    if remove_temp_outputs:
        os.system('rm cmd_modeller_ini.py')
    
    # Calculating all Atom contacts
    if 'contact-chainID_allAtoms' not in os.listdir('../../../modelling_scripts/'):
        os.system('g++ ../../../modelling_scripts/contact-chainID_allAtoms.cpp -o ../../modelling_scripts/contact-chainID_allAtoms')
    os.popen('../../../modelling_scripts/contact-chainID_allAtoms %s.ini %s > all_contacts_%s.list' %(target_name.upper(), cutoff, template_ID)).read()
    
    ### Writing anchors contact list ###
    
    with open( 'all_contacts_%s.list' %template_ID, 'r') as contacts:
        with open('contacts_%s.list' %template_ID, 'w') as output:
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
    
    if remove_temp_outputs:
        os.system('rm all_contacts_%s.list' %template_ID)
    
    with open('MyLoop.py', 'w') as myloopscript:
        MyL_temp = open('../../../modelling_scripts/MyLoop_template.py', 'r')
        for line in MyL_temp:
            if 'self.residue_range' in line:
                myloopscript.write(line %(anchor_1 + 2, anchor_2))
            elif 'contact_file = open' in line:
                myloopscript.write(line %template_ID)
            else:
                myloopscript.write(line)
        MyL_temp.close()
        
    #    with open('instructions.txt', 'w') as instr_file:
    #        instr_file.write(template_ID + ' ' + str(anch_1) + ' ' + str(anch_2) + ' ' + str(modeller_renum))
    
    #Finally launching Modeller. Hopefully.
    
    t1 = time.time()
    
    with open('cmd_modeller.py', 'w') as modscript:
        cmd_m_temp = open('../../../modelling_scripts/cmd_modeller_template.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line %final_alifile_name)
            elif 'knowns' in line:
                modscript.write(line %(template_ID, target_name.upper()))
            else:
                modscript.write(line)
        cmd_m_temp.close()
    
    os.popen('python3 cmd_modeller.py > modeller.log').read()
    
    '''
    with open('instructions.txt', 'w') as instr_file:
        instr_file.write(template_ID + ' ' + str(anchor_1) + ' ' + str(anchor_2))
    
    t1 = time.time()
    
    os.popen('python2.7 ../../modelling_scripts/cmd_modeller.py %s %s %s' %(final_alifile_name, template_ID, target_name.upper())).read()
    '''
    
    t2 = time.time()
    tf = t2 - t1
    
    print('The modelling took %i seconds' %tf)
    
    os.popen('python3 ../../../modelling_scripts/get_dope_scores.py modeller.log DOPE.tsv %s %s' %(template_ID, target_name))

