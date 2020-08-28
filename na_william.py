#TODO:
# - Prepare template MHC seqs?
# - Define Mutation case
# - Download fasta seq
# - Run NetMHCpan to identify putative binding peptides from mutated
# - Use peptide seq to identify one template from each allele
# - Do the models (12 + 11 + 10 + 9 + 8 = 50)

# Usage: python na_william.py <outdir name> <taskID> <n_tasks> <n_cores>

from modelling_scripts import data_prep
import pickle
from Bio import SeqIO
#import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
#from Bio.PDB import PDBParser
import os
import time
from Bio.SubsMat import MatrixInfo
from random import choice
import csv
import sys
from joblib import Parallel, delayed
import copy


with open ('./william/mutations.txt', 'r') as mutfile:
    proteins = []
    mutations = []
    for line in mutfile:
        protein = line.split(' ')[0]
        mutation = line.split(' ')[1]
        proteins.append(protein)
        mutations.append((protein, mutation))
proteins = list(set(proteins))
'''
os.chdir('./william/fasta')
for protein in proteins:
    os.system('wget https://www.uniprot.org/uniprot/%s.fasta' %protein)
os.chdir('../../')
'''

###########################################
### Organizing dicts ###

IDD_file = open('data/csv_pkl_files/IDs_ChainsCounts_dict.pkl', 'rb')
IDD = pickle.load(IDD_file)
bad_IDs = pickle.load(IDD_file)
IDD_file.close()

### Organizing Allele IDs in a dictionary ###
allele_ID = {}
for key in IDD:
    #if 'HLA' in IDD[key]['allele']:
    try:
        allele_ID[IDD[key]['allele']].append(key)
    except KeyError:
        allele_ID[IDD[key]['allele']] = [key]

### Retrieving common HLA alleles between templates and netMHCpanI allele list ###
with open('/home/dariomarzella/NetMHCpanI/netMHCpan-4.1/Linux_x86_64/data/MHC_pseudo.dat') as infile:
    NMP_allele_set = []
    for line in infile:
        NMP_allele_set.append(line.split(' ')[0])

common_alleles = []
for allele in allele_ID:
    al = allele.replace('*', '')
    if al.startswith('HLA-') and al in NMP_allele_set:
        common_alleles.append((allele, al))
print(common_alleles, len(common_alleles))

leads = []
for case in os.listdir('./william/fasta/targets/'):
    NMP_outdir = case.split('.')[0]
    try:
        os.mkdir('./william/fasta/netMHCpanI_outputs/' + NMP_outdir)
    except FileExistsError:
        pass


    for allele in common_alleles:
        os.system('$netMHCpanI -a %s ./william/fasta/targets/%s > ./william/fasta/netMHCpanI_outputs/%s/%s_%s'
                  %(allele[1], case, NMP_outdir, allele[1], case)) #allele[1].replace('-','').replace(':','')))
        print(allele[1] + ' DONE')

    hits = []
    for filename in os.listdir('./william/fasta/netMHCpanI_outputs/' + NMP_outdir):
        with open('./william/fasta/netMHCpanI_outputs/' + NMP_outdir + '/' + filename) as infile:
            for line in infile:
                if '<=' in line:
                    fline = [x for x in line.split(' ') if x != '']
                    hits.append((filename, fline))
                    if int(fline[0]) < 12 and (int(fline[0]) + len(fline[9])) > 12:
                        #print(filename, fline)
                        #print('MUTATION IN PEPTIDE')
                        leads.append((fline[9], fline[1], NMP_outdir))

    print('case %s LEADS\n' %case)
    print(leads)
######################################################################################################
######################################################################################################
######################################################################################################


start_time = time.time()
###########################################
### Organizing dicts ###

###########################################
###  Setting variables

outdir_name = sys.argv[1]

'''
if outdir_name == 'wildtype':
    pept_seqs = data_prep.get_peptides_from_csv('data/csv_pkl_files/table_1_AnAnalysisofNaturalTCellResponsestoPredictedTumorNeoepitopes.csv', 3, 4, ',')
elif outdir_name == 'neoantigens':
    pept_seqs = data_prep.get_peptides_from_csv('data/csv_pkl_files/table_1_AnAnalysisofNaturalTCellResponsestoPredictedTumorNeoepitopes.csv', 1, 4, ',')
else:
    raise Exception('Invalid output directory name')
'''
pept_seqs = leads


###
#pept_seqs = pept_seqs[:4]
###

'''
try:
    pepts_start = int(sys.argv[2])
except:
    pepts_start = 0
try:
    pepts_end = int(sys.argv[3])
except:
    pepts_end = len(pept_seqs)
'''
taskID = int(sys.argv[2])
n_tasks = int(sys.argv[3])
n_cores = int(sys.argv[4])

step_width = len(pept_seqs)/n_tasks
pepts_start = int(step_width * (taskID-1))
pepts_end = int(step_width * taskID)

print('######################################')
print(taskID)
print(n_tasks)
print(pepts_start)
print(pepts_end)
print('######################################')

anch_dict = { 8: (1, 7), 9 : (1, 8), 10 : (1, 9),
              11 : (1, 10), 12 : (1, 11)}

remove_temp_outputs = False
cutoff = 5

maxs = []
maxsl = []
non_modelled = []
print_results = False

filename_start = 'BL00'
filename_end = '0001.pdb'

###########################################

def na_model(k, pept_seq, outdirname_add = None):
    #for k, pept_seq in zip(range(pepts_start, pepts_end), pept_seqs[pepts_start:pepts_end]):
    t1 = time.time()
    ext_flag = False
    query = 'query_' + str(k + 1)
    print( '## Wokring on %s ##' %query)
    '''
    for folder in os.listdir('outputs/%s' %outdir_name):
        if 'query' in folder:
            if str(k+1) == folder.split('_')[2]:
                print('WARNING: Existing directory for this query.')
                for name in os.listdir('outputs/%s/%s' %(outdir_name, folder)):
                    if filename_start in name and filename_end in name:
                        print('WARNING: This query has already been modelled. Moving on.')
                        ext_flag = True
                        break
                if ext_flag == True:
                    break
    if ext_flag == True:
        print('Exiting here.')
        return None
    '''
    pept = pept_seq[0]
    allele = pept_seq[1]
    outdirname_add = pept_seq[2]
    #outdirname_add = copy.deepcopy(allele)
    length = len(pept)

    max_pos = -1000
    pos_list = []

    anch_1, anch_2 = anch_dict[length]

    homolog_allele = '--NONE--'
    if allele.startswith('HLA'):      # Human
        if allele in allele_ID.keys():
            pass
        elif allele[:9] in allele_ID.keys():
            allele = allele[:9]
        else:
            allele = allele[:6]
    elif allele.startswith('H2'):    # Mouse
        #homolog_allele = 'RT1'
        if allele in allele_ID.keys():
            pass
        elif allele[:4] in allele_ID.keys():
            allele = allele[:4]
        else:
            allele = allele[:3]
    elif allele.startswith('RT1'):          # Rat
        homolog_allele = 'H2'
        if allele in allele_ID.keys():
            pass
        elif allele[:5] in allele_ID.keys():
            allele = allele[:5]
        else:
            allele = allele[:4]
    elif allele.startswith('BoLA'):        # Bovine
        if allele in allele_ID.keys():
            pass
        elif allele[:10] in allele_ID.keys():
            allele = allele[:10]
        elif allele[:7] in allele_ID.keys():
            allele = allele[:7]
        else:
            allele = allele[:5]
    elif allele.startswith('SLA'):        # Suine
        if allele in allele_ID.keys():
            pass
        elif allele[:9] in allele_ID.keys():
            allele = allele[:9]
        elif allele[:6] in allele_ID.keys():
            allele = allele[:6]
        else:
            allele = allele[:4]
    elif allele.startswith('BF2'):        # Chicken
        if allele in allele_ID.keys():
            pass
        elif allele[:6] in allele_ID.keys():
            allele = allele[:6]
        else:
            allele = allele[:4]
    elif allele.startswith('Mamu'):       # Monkey
        if allele in allele_ID.keys():
            pass
        elif allele[:13] in allele_ID.keys():
            allele = allele[:13]
        elif allele[:9] in allele_ID.keys():
            allele = allele[:9]
        else:
            allele = allele[:5]
    elif allele.startswith('Eqca'):        # Horse
        if allele in allele_ID.keys():
            pass
        elif allele[:10] in allele_ID.keys():
            allele = allele[:10]
        elif allele[:7] in allele_ID.keys():
            allele = allele[:7]
        else:
            allele = allele[:5]

    for ID in IDD:
        if allele in IDD[ID]['allele'] or homolog_allele in IDD[ID]['allele']:                       ## Same Allele
            score = 0
            temp_pept = IDD[ID]['pept_seq']
            min_len = min([length, len(temp_pept)])
            #score -= (abs(length - len(temp_pept)) * 20)      ## Penalty for gaps
            score -= int(2.4 ** (2 + abs(length - len(temp_pept))))
            for (aa, bb) in zip(pept[:min_len], temp_pept[:min_len]):
                try:
                    score += MatrixInfo.pam30[aa, bb]
                except KeyError:
                    try:
                        score += MatrixInfo.pam30[bb, aa]
                    except KeyError:
                        #print('Broken peptide in IDD. Passing over.')
                        score = -50
                        pass

            if score > max_pos:
                max_pos = score
            pos_list.append((score, temp_pept, ID))
        else:
            pass

    max_list = []
    for pos in pos_list:
        if pos[0] == max_pos:
            max_list.append(pos)
    maxs.append(max_pos)
    maxsl.append(len(max_list))

    if len(max_list) == 0:
        return (query, pept, allele, "NA", "NA", "NA", "No positive scoring template peptides")
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

    ###################################
    #   CHANGING WORKING DIRECTORY
    ###################################

    sequences, empty_seqs = data_prep.get_pdb_seq([template_ID])

    outdir = ('outputs/%s/%s_%s_%s_%s' %(outdir_name, outdirname_add, allele, template_ID.lower(), query))
    try:
        os.mkdir(outdir)
    except FileExistsError:
        print('WARNING: Existing directory.')
        #continue
    os.chdir(outdir)

    #Preparing template and target sequences for the .ali file, launching Muscle

    template_seqr = SeqRecord(Seq(sequences[0]['M'], IUPAC.protein), id=template_ID, name = template_ID)
    target_seqr = SeqRecord(Seq(sequences[0]['M'], IUPAC.protein), id=query, name = query)
    SeqIO.write((template_seqr, target_seqr), "%s.fasta" %template_ID, "fasta")

    os.system('muscle -in %s.fasta -out %s.afa -quiet' %(template_ID, template_ID))
    os.system('rm %s.fasta' %template_ID)

    ali_template_pept = copy.deepcopy(template_pept)
    ali_target_pept = copy.deepcopy(pept)
    if len(template_pept) > length:
        ali_target_pept = pept[0:5] + ('-'*(len(template_pept)-length)) + pept[5:]
    elif length > len(template_pept):
        ali_template_pept = template_pept[0:5] + ('-'*(length-len(template_pept))) + template_pept[5:]
    final_alifile_name = '%s.ali' %template_ID
    final_alifile = open(final_alifile_name, 'w')
    i = 0
    target_ini_flag = False
    template_ini_flag = False
    modeller_renum = 1
    for line in open('%s.afa' %template_ID, 'r'):
        #print(line)
        if line.startswith('>') and i == 0:
            final_alifile.write('>P1;' + line.split(' ')[0].strip('>') + '\n')
            final_alifile.write('structure:../../../data/PDBs/%s_MP.pdb:%s:M:%s:P::::\n' %(template_ID, str(sequences[0]['M_st_ID']), str(len(ali_template_pept))))
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

    #with open('instructions.txt', 'w') as instr_file:
    #    instr_file.write(template_ID + ' ' + str(1) + ' ' + str(9))

    '''
    os.system('pdb_splitchain ../../data/PDBs/%s_MP.pdb' %template_ID)
    for chain in ['M', 'P']:
        os.system('mv %s_MP_%s.pdb ../../data/PDBs/%s_MP_%s.pdb' %(template_ID, chain, template_ID, chain))
        os.system('pdb_reres -1 ../../data/PDBs/%s_MP_%s.pdb > ../../data/PDBs/%s_MP_%s.pdb.renum' %(template_ID, chain, template_ID, chain))
    os.system('pdb_merge ../../data/PDBs/%s_MP_*.pdb.renum > ../../data/PDBs/%s_MP_reres.pdb' %(template_ID, template_ID))
    os.system('rm -f data/PDBs/%s_MP_[A-Z].pdb' %template_ID)
    os.system('rm -f ../../data/PDBs/*.pdb.renum')
    '''

    with open('MyLoop.py', 'w') as myloopscript:
        MyL_temp = open('../../../modelling_scripts/MyLoop_template.py', 'r')
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
        cmd_m_temp = open('../../../modelling_scripts/cmd_modeller_ini.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line %final_alifile_name)
            elif 'knowns' in line:
                modscript.write(line %(template_ID, query))
            else:
                modscript.write(line)
        cmd_m_temp.close()

    os.popen('python3 cmd_modeller_ini.py').read()

    if remove_temp_outputs:
        os.system('rm cmd_modeller_ini.py')

    # Calculating all Atom contacts
    if "contact-chainID_allAtoms" not in os.listdir('../../../modelling_scripts'):
        os.popen('g++ ../../../modelling_scripts/contact-chainID_allAtoms.cpp -o ../../../modelling_scripts/contact-chainID_allAtoms').read()
    os.popen('../../../modelling_scripts/contact-chainID_allAtoms %s.ini %s > all_contacts_%s.list' %(query, cutoff, template_ID)).read()

    #Selecting only the anchors contacts

    real_anchor_2 = None
    anch_1_same = False
    anch_2_same = False

    if template_pept[anch_1] == pept[anch_1]:
        anch_1_same = True
    if len(template_pept) >= length:
        if template_pept[anch_2] == pept[anch_2]:
            anch_2_same = True
    #else:
    #    if template_pept[-1] == pept[anch_2]:
    #        anch_2_same = True

    if (anch_2 + 1) == len(template_pept):
        pass
    elif (anch_2 + 1) > len(template_pept):
        real_anchor_2 = len(template_pept)  #????
    elif (anch_2 + 1) < len(template_pept):
        real_anchor_2 = len(template_pept)

    if pept[anch_1] == 'G':
        first_gly = True
    else:
        first_gly= False
    if pept[-1] == 'G':
        last_gly = True
    else:
        last_gly = False

    ### Writing anchors contact list ###

    with open( 'all_contacts_%s.list' %template_ID, 'r') as contacts:                             # Template contacts
        with open('contacts_%s.list' %template_ID, 'w') as output:
            if real_anchor_2:                                                                     ### If the target peptide is longer than the templtate peptide #TODO: This can be removed?
                for line in contacts:
                    #print(line[30:33])
                    p_aa_id = line.split("\t")[7]                                                 ### position id of the template peptide residue
                    p_atom = line.split("\t")[8]                                                  ### atom name of the template peptide residue
                    #m_aa_id = (line.split("\t")[2]).split(' ')[0]
                    if anch_1_same == True:                                                       ### If the target anchor 1 residue is the same as the template anchor 1 residue
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
                    '''
                    else:
                        if int(p_aa_id) == (anch_1+1) and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line)
                    '''
                    if anch_2_same == True:                                                       ### If the target anchor 2 residue is the same as the template anchor 2 residue
                        if int(p_aa_id) == real_anchor_2:
                            output.write(line[:30] + str(anch_2+1) + line[34:])
                    else:
                        if int(p_aa_id) == (anch_2+1):
                            if last_gly:
                                if 'CA' in p_atom:
                                    output.write(line[:30] + str(anch_2+1) + line[34:])
                            else:
                                if 'CA' in p_atom or 'CB' in p_atom:
                                    output.write(line[:30] + str(anch_2+1) + line[34:])
                    '''
                    else:
                        if int(p_aa_id) == real_anchor_2 and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line[:30] + str(anch_2+1) + line[34:])
                    '''
            else:
                for line in contacts:
                    #print(line.split("\t"))
                    p_aa_id = line.split("\t")[7]
                    p_atom = line.split("\t")[8]
                    #m_aa_id = (line.split("\t")[2]).split(' ')[0]
                    if anch_1_same == True:                                                       ### If the target anchor 1 residue is the same as the template anchor 1 residue
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
                    '''
                    else:
                        if int(p_aa_id) == (anch_1+1) and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line)
                    '''
                    if anch_2_same == True:                                                       ### If the target anchor 2 residue is the same as the template anchor 2 residue
                        if int(p_aa_id) == real_anchor_2:
                            output.write(line[:30] + str(anch_2+1) + line[34:])
                    else:
                        if int(p_aa_id) == (anch_2+1):
                            if last_gly:
                                if 'CA' in p_atom:
                                    output.write(line)
                            else:
                                if 'CA' in p_atom or 'CB' in p_atom:
                                    output.write(line)

    if remove_temp_outputs:
        os.system('rm all_contacts_%s.list' %template_ID)

    with open('MyLoop.py', 'w') as myloopscript:
        MyL_temp = open('../../../modelling_scripts/MyLoop_template.py', 'r')
        for line in MyL_temp:
            if 'self.residue_range' in line:
                myloopscript.write(line %(anch_1 + 2, anch_2))
            elif 'contact_file = open' in line:
                myloopscript.write(line %template_ID)
            else:
                myloopscript.write(line)
        MyL_temp.close()

    #    with open('instructions.txt', 'w') as instr_file:
    #        instr_file.write(template_ID + ' ' + str(anch_1) + ' ' + str(anch_2) + ' ' + str(modeller_renum))

    #Finally launching Modeller. Hopefully.

    smt = time.time()

    with open('cmd_modeller.py', 'w') as modscript:
        cmd_m_temp = open('../../../modelling_scripts/cmd_modeller_template.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line %final_alifile_name)
            elif 'knowns' in line:
                modscript.write(line %(template_ID, query))
            else:
                modscript.write(line)
        cmd_m_temp.close()

    os.popen('python3 cmd_modeller.py > modeller.log').read()

    t2 = time.time()
    mt = t2 - smt
    tf = t2 - t1

    print('The modelling took %i seconds' %mt)
    print('The whole case took %i seconds' %tf)

    if mt < 3:
        os.chdir('../../../')
        return (query, pept, allele, template_ID, template_pept, IDD[template_ID]['allele'], 'Something gone wrong in the modelling')

    os.chdir('../../../')
    return None

######################################################################################################
######################################################################################################
non_modelled = Parallel(n_jobs = n_cores)(delayed(na_model)(k, pept_seq) for k, pept_seq in enumerate(pept_seqs[pepts_start:pepts_end]))
non_modelled = filter(None, non_modelled)

with open('outputs/%s/non_modelled_%s.csv' %(outdir_name, taskID), 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['QUERY', 'NEOANTIGEN', 'QUERY ALLELE', 'TEMPLATE ID', 'TEMPLATE PEPTIDE', 'TEMPLATE ALLELE', 'ERROR'])
    for query in non_modelled:
        tsv_writer.writerow(query)

final_time = time.time()
total_time = final_time - start_time
print('The whole pipeline took: ', total_time)
