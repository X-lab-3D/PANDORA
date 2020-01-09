
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

#IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/final_mhc1_3d_structure_data_with_pdb_ids.tsv')

## TODO:
## Retrieve IDs : Alleles
## Read user's input (seq and allele)
## Select template structure and get its sequence
## Compare User's input with same allele sequences (select with identity the best one)

### Retriving Dictionary with PDB IDs and chain lengths ###
IDD_file = open('data/IDs_ChainsCounts_dict.pkl', 'rb')
IDD = pickle.load(IDD_file)
bad_IDs = pickle.load(IDD_file)
IDD_file.close()

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
print('Please select your sequence file.')
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
sequences = data_prep.get_pdb_seq(allele_ID[allele])
ID_seqs_dict = {}
for i, ID in enumerate(allele_ID[allele]):
    template = SeqRecord(Seq(sequences[i]['M'], IUPAC.protein), id=ID, name = allele + ID)
    SeqIO.write((template, A_target), "data/FASTAs/%s.fasta" %ID, "fasta")
    ID_seqs_dict[ID] = ((sequences[i]['P'], len(sequences[i]['M'])))

### Identifying the template ID and its peptide sequence ### 
if len(allele_ID[allele]) == 1:   ### In case we have only one template for this allele ###
    template_ID = allele_ID[allele][0]
    peptide_seq = ID_seqs_dict[template_ID][0]
    #subprocess.check_call(['muscle', '-in %s.fasta -clw' %template_ID])
    t1 = time.time()
    os.system('muscle -in data/FASTAs/%s.fasta -out data/Alignments/%s.afa' %(ID, ID))
    tf = time.time() - t1
    print('IT TOOK: ')
    print(tf)
else:                             ### In case we have multiple templates for this allele ###
    score_dict = {}
    max_id = 0
    putative_templates = []
    for ID in allele_ID[allele]:
        alifile = 'data/Alignments/%s.afa' %ID
        muscle_commands = ['muscle', '-in', 'data/FASTAs/%s.fasta' %ID, '-out', '%s' %alifile, '-quiet']
        #os.system('muscle -in data/FASTAs/%s.fasta -out %s' %(ID, alifile))
        subprocess.check_call(muscle_commands)
        os.system('rm data/FASTAs/%s.fasta' %ID)
        aligned = list(SeqIO.parse(alifile, 'fasta'))
        id_count = 0
        for x, y in zip((aligned[0].seq), str(aligned[1].seq)):
            if x == y:
                id_count += 1
            else:
                pass
        if id_count > max_id:
            max_id = id_count
        score_dict[ID] = id_count
        #print('ID COUNT:  ', id_count)
    for key in score_dict:
        if score_dict[key] == max_id:
            putative_templates.append(key)
        else:
            os.system('rm data/Alignments/%s.afa' %key)
    templates_dict = {}
    for i, template in enumerate(putative_templates):
        templates_dict[i] = template
    print('##############################')
    print('Please select one template between the ones listed below')
    print('##############################')
    print('')
    print('')
    print(templates_dict)
    template_ID = templates_dict[int(input())]
    #template_ID = templates_dict[3]
    #template_ID = choice(putative_templates)
    peptide_seq = ID_seqs_dict[template_ID][0]
    
    print('##############################')
    print('Please input one restrain distance cutoff')
    print('##############################')
    print('')
    cutoff = 5

### Writing a final .ali file with Template sequence / template pept sequence ; Target sequence / target pept sequence ###
    
os.system('pdb_reres -1 data/PDBs/%s_MP.pdb > data/PDBs/%s_MP_reres.pdb' %(template_ID, template_ID))  # Renumbering the residues

final_alifile_name = 'data/Alignments/%s.ali' %template_ID
final_alifile = open(final_alifile_name, 'w')
i = 0
for line in open('data/Alignments/%s.afa' %template_ID, 'r'):
    #print(line)
    if line.startswith('>') and i == 0:
        final_alifile.write('>P1;' + line.split(' ')[0].strip('>') + '\n')
        final_alifile.write('structure:data/PDBs/%s_MP_reres.pdb:1:M:9:P::::\n' %template_ID)
        i += 1
    elif line.startswith('>') and i == 1:
        final_alifile.write('/' + peptide_seq + '*')
        final_alifile.write('\n')
        final_alifile.write('\n>P1;' + line.split(':')[0].strip('>') + '\n')
        final_alifile.write('sequence:::::::::\n')
    else:
        final_alifile.write(line.rstrip())
final_alifile.write('/' + str(P_target.seq) + '*')
final_alifile.close()

#final_alifile_name = 'data/Alignments/1k5n_1ogt.ali'

os.system('rm data/Alignments/*.afa')

# Calculating all Atom contacts
os.popen('modelling_scripts/contact-chainID_allAtoms data/PDBs/%s_MP_reres.pdb %s > data/all_contacts_%s.list' %(template_ID, cutoff, template_ID)).read()

#Selecting only the anchors contacts
print('##############################')
print('Please input, one per time, the anchor positions (only the int number)')
print('##############################')
print('')
print('The peptide of the selected template is %i residue long. The target one is %i residue long' %(len(str(peptide_seq)), len(str(P_target.seq))))
print('')
print('The template peptide sequence is %s, the target peptide seuqnce is %s' %(str(peptide_seq), str(P_target.seq)))
print('ANCHOR 1:')
anchor_1 = int(input()) + ID_seqs_dict[template_ID][1]
print('ANCHOR 2:')
anchor_2 = int(input()) + ID_seqs_dict[template_ID][1]
real_anchor_2 = None

if anchor_2 > (len(str(peptide_seq)) + ID_seqs_dict[template_ID][1]):
    real_anchor_2 = (len(str(peptide_seq)) + ID_seqs_dict[template_ID][1])
'''
print('')
print('############################################################')
print('WARNING: you are selecting all the restrains, not only CA or CB')
print('############################################################')
print('')
'''
with open( 'data/all_contacts_%s.list' %template_ID, 'r') as contacts:
    with open('data/contacts_P%i_P%i.list' %(anchor_1, anchor_2), 'w') as output:
        if real_anchor_2:
            for line in contacts:
                #print(line[30:33])
                p_aa_id = line.split("\t")[7]
                p_atom = line.split("\t")[8]
                m_aa_id = (line.split("\t")[2]).split(' ')[0]
                if int(p_aa_id) == anchor_1 and ('CA' in p_atom or 'CB' in p_atom):
                    output.write(line)
                elif int(p_aa_id) == real_anchor_2 and ('CA' in p_atom or 'CB' in p_atom):
                    output.write(line[:30] + str(anchor_2) + line[34:])
        else:
            for line in contacts:
                #print(line.split("\t"))
                p_aa_id = line.split("\t")[7]
                p_atom = line.split("\t")[8]
                m_aa_id = (line.split("\t")[2]).split(' ')[0]
                if (int(p_aa_id) == anchor_1 or int(p_aa_id) == anchor_2) and ('CA' in p_atom or 'CB' in p_atom):
                    output.write(line)

            
with open('data/instructions.txt', 'w') as instr_file:
    instr_file.write(str(anchor_1) + '\n')
    instr_file.write(str(anchor_2))
    
#Finally launching Modeller. Hopefully.
                
#command = ['python', 'modelling_scripts/cmd_modeller.py', final_alifile_name, template_ID, '>3ROO:A']

#proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
#(out, err) = proc.communicate()
                
#raise Exception('OK.')
os.popen('/usr/bin/python2.7 modelling_scripts/cmd_modeller.py %s %s %s' %(final_alifile_name, template_ID, target_name.upper())).read()

#os.popen('/usr/bin/python2.7 modelling_scripts/cmd_modeller.py %s 1k5n_MP query_1ogt_MP' %final_alifile_name).read()





