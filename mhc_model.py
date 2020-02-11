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

###IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/final_mhc1_3d_structure_data_with_pdb_ids.tsv')

#IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/mhcI_structures_IEDB.csv', 42, 43, ';', [0,1])

### Retriving Dictionary with PDB IDs and chain lengths ###

IDD_file = open('data/csv_pkl_files/IDs_ChainsCounts_dict.pkl', 'rb')
IDD = pickle.load(IDD_file)
bad_IDs = pickle.load(IDD_file)
IDD_file.close()

remove_temp_outputs = True
cutoff = 5

### Organizing Allele IDs in a dictionary ###
allele_ID = {}
for key in IDD:
    try:
        allele_ID[IDD[key]['allele']].append(key)
    except KeyError:
        allele_ID[IDD[key]['allele']] = [key]

### Retriving target allele and sequence file ###

print('##############################')
print('Please enter your allele name. You can choose between the ones listed below')
print('##############################')
print('')
print('')
print(allele_ID.keys())
allele = input()
#allele, inseq_file = 'FLA-E*01801', 'data/5xmf.fasta'
#allele, inseq_file = 'HLA-B*27:09', 'data/1ogt.fasta'

print('##############################')
print('Please select your target sequence file.')
print('##############################')
print('')

for file in glob.glob(os.getcwd() + '/data/Targets/*.fasta'):
    print(file.split('/')[-1].split('.')[0])

#inseq_file = 'data/5xmf.fasta'
#inseq_file = 'data/1tom.fasta'
target_name = input()
inseq_file = 'data/Targets/%s.fasta' %target_name
inseqs = list(SeqIO.parse(inseq_file, 'fasta'))
A_target = inseqs[0]
P_target = inseqs[1]

### Producing a Fasta file for each putetive template-target alignment ###
sequences, empty_seqs = data_prep.get_pdb_seq(allele_ID[allele])
for bad_ID in empty_seqs:
    del IDD[bad_ID]
    allele_ID[allele].remove(bad_ID)

ID_seqs_dict = {}

### Identifying the template ID and its peptide sequence ###
if len(allele_ID[allele]) == 1:   ### In case we have only one template for this allele ###
    template_ID = allele_ID[allele][0]
    outdir = ('outputs/%s_%s' %(template_ID.lower(), target_name))
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
    template_pept = ID_seqs_dict[template_ID][0]
    #subprocess.check_call(['muscle', '-in %s.fasta -clw' %template_ID])
    t1 = time.time()
    os.system('muscle -in %s.fasta -out %s.afa' %(template_ID, template_ID))
    os.system('rm %s.fasta' %template_ID)
    tf = time.time() - t1
    print('IT TOOK: ')
    print(tf)
else:                             ### In case we have multiple templates for this allele ###
    score_dict = {}
    max_id = 0
    templates_dict = allele_ID[allele]
    print('##############################')
    print('Please select one template between the ones listed below. You can enter the key number or the ID, if you are sure it is present')
    print('##############################')
    print('')
    print('')
    for i, j in enumerate(allele_ID[allele]):
        print(i,j)
    temp_in = input()
    try:
        template_ID = allele_ID[allele][int(temp_in)]
    except:
        template_ID = temp_in
    outdir = ('outputs/%s_%s' %(template_ID.lower(), target_name))
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
    os.system('rm %s.fasta' %template_ID)

    template_pept = ID_seqs_dict[template_ID][0]

### Writing a final .ali file with Template sequence / template pept sequence ; Target sequence / target pept sequence ###

#os.system('pdb_reres -1 data/PDBs/%s_MP.pdb > data/PDBs/%s_MP_reres.pdb' %(template_ID, template_ID))  # Renumbering the residues

'''
os.system('pdb_splitchain data/PDBs/%s_MP.pdb' %template_ID)
for chain in ['M', 'P']:
    os.system('mv %s_MP_%s.pdb data/PDBs/%s_MP_%s.pdb' %(template_ID, chain, template_ID, chain))
    os.system('pdb_reres -1 data/PDBs/%s_MP_%s.pdb > data/PDBs/%s_MP_%s.pdb.renum' %(template_ID, chain, template_ID, chain))
os.system('pdb_merge data/PDBs/%s_MP_*.pdb.renum > data/PDBs/%s_MP_reres.pdb' %(template_ID, template_ID))
os.system('rm -f data/PDBs/%s_MP_[A-Z].pdb' %template_ID)
os.system('rm data/PDBs/*.pdb.renum')
'''

final_alifile_name = '%s.ali' %template_ID
final_alifile = open(final_alifile_name, 'w')
i = 0
target_ini_flag = False
template_ini_flag = False

for line in open('%s.afa' %template_ID, 'r'):
    #print(line)
    if line.startswith('>') and i == 0:
        final_alifile.write('>P1;' + line.split(' ')[0].strip('>') + '\n')
        final_alifile.write('structure:../../data/PDBs/%s_MP.pdb:1:M:9:P::::\n' %template_ID)
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

with open('instructions.txt', 'w') as instr_file:
    instr_file.write(template_ID + ' ' + str(1) + ' ' + str(9))
os.popen('/usr/bin/python2.7 ../../modelling_scripts/cmd_modeller_ini.py %s %s %s' %(final_alifile_name, template_ID, target_name.upper())).read()

# Calculating all Atom contacts
os.popen('../../modelling_scripts/contact-chainID_allAtoms %s.ini %s > all_contacts_%s.list' %(target_name.upper(), cutoff, template_ID)).read()

#Selecting only the anchors contacts
print('##############################')
print('Please input, one per time, the anchor positions (only the int number)')
print('##############################')
print('')
print('The peptide of the selected template is %i residue long. The target one is %i residue long' %(len(str(template_pept)), len(str(P_target.seq))))
print('')
print('The template peptide sequence is %s, the target peptide seuqnce is %s' %(str(template_pept), str(P_target.seq)))
print('ANCHOR 1:')
anchor_1 = int(input())
#anchor_1 = int(anch_1_in) + ID_seqs_dict[template_ID][1]
print('ANCHOR 2:')
anchor_2 = int(input())
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

t1 = time.time()
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

with open('instructions.txt', 'w') as instr_file:
    instr_file.write(template_ID + ' ' + str(anchor_1) + ' ' + str(anchor_2))
    
os.popen('/usr/bin/python2.7 ../../modelling_scripts/cmd_modeller.py %s %s %s' %(final_alifile_name, template_ID, target_name.upper())).read()


t2 = time.time()
tf = t2 - t1

print('It took %i seconds' %tf)

