
#!/usr/bin/python
###    ###
from modelling_scripts import data_prep
import pickle
from Bio import SeqIO
#import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
import time
from Bio.SubsMat import MatrixInfo
from random import choice
import csv
#import sys

### Retriving Dictionary with PDB IDs and chain lengths ###

##IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/final_mhc1_3d_structure_data_with_pdb_ids.tsv')
#IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/csv_pkl_files/mhcI_structures_IEDB.csv', 42, 43, ';', empty_rows=[0,1])

#outdir_name = sys.argv[1]

start_time = time.time()
outdir_name = 'benchmark'

###########################################
### Organizing dicts ###

IDD_file = open('data/csv_pkl_files/IDs_ChainsCounts_dict.pkl', 'rb')
main_IDD = pickle.load(IDD_file)
bad_IDs = pickle.load(IDD_file)
IDD_file.close()

raise Exception('OK')

'''
### Organizing Allele IDs in a dictionary ###
allele_ID = {}
for key in IDD:
    if 'HLA' in IDD[key]['allele']:
        try:
            allele_ID[IDD[key]['allele']].append(key)
        except KeyError:
            allele_ID[IDD[key]['allele']] = [key]
'''
###########################################
###  Setting variables 
            
# tsv with: Target ID, pept seq, allele
# pept_seqs = data_prep.get_peptides_from_csv('benchmark/cross_validation_set.tsv', 3, 4, '\t') ### This can be the input requested for the benchmark function


pept_seqs = []
with open('benchmark/PANDORA_benchmark_dataset.tsv', 'r') as peptsfile:
    spamreader = csv.reader(peptsfile, delimiter='\t')
    for i, row in enumerate(spamreader):
        if i == 0:
            pass
        else:
            target_id = row[0]
            seq = row[1]
            allele = row[2]
            if 'HLA' in allele:
                star_allele = (allele[0:5]+'*'+allele[5:])
                pept_seqs.append((target_id, seq, star_allele))
            else:
                pept_seqs.append((target_id, seq, allele))

anch_dict = { 8: (1, 7), 9 : (1, 8), 10 : (1, 9),   #Python numbering 0-X
              11 : (1, 10), 12 : (1, 11)}

remove_temp_outputs = False
cutoff = 5

maxs = []
maxsl = []
non_modelled = []
print_results = False

filename_start = 'BL00'
filename_end = '0001.pdb'
best_rmsds = {}

###########################################

for k, pept_seq in enumerate(pept_seqs):
    t1 = time.time()
    
    target_id = pept_seq[0]
    pept = pept_seq[1]
    allele = pept_seq[2]
    length = len(pept)
    M_chain = data_prep.get_seqs('./data/PDBs/%s_MP.pdb' %target_id)['M']
    
    if length < 8 or length > 12:
        print('###')
        print('Invalid Peptide Length. Exiting here.')
        print('###')
        continue
    
    IDD = main_IDD
    del IDD[target_id]
    
    ### Organizing Allele IDs in a dictionary ###
    allele_ID = {}
    for key in IDD:
        if 'HLA' in IDD[key]['allele']:
            try:
                allele_ID[IDD[key]['allele']].append(key)
            except KeyError:
                allele_ID[IDD[key]['allele']] = [key]
    
    ext_flag = False
    query = 'query_' + str(k + 1)
    print( '## Wokring on %s, structure %s ##' %(query, target_id))
    for folder in os.listdir('outputs/%s' %outdir_name):
        if target_id in folder:
            print('WARNING: Existing directory for this query.')
            ext_flag = True
            break
    if ext_flag == True:
        print('Exiting here. This is only a DEBUG message')
        continue
        
    max_pos = -200
    pos_list = []
    
    anch_1, anch_2 = anch_dict[length]
    
    if allele.startswith('HLA'):      # Human
        if allele in allele_ID.keys():
            pass
        elif allele[:9] in allele_ID.keys():
            allele = allele[:9]
        else:
            allele = allele[:6]
    elif allele.startswith('H2'):    # Mouse
        if allele in allele_ID.keys():
            pass
        elif allele[:4] in allele_ID.keys(): 
            allele = allele[:4]
        else:
            allele = allele[:3]
    elif allele.startswith('RT1'):
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
        if allele in IDD[ID]['allele']:                       ## Same Allele
            score = 0
            temp_pept = IDD[ID]['pept_seq']
            min_len = min([length, len(temp_pept)])
            score -= (abs(length - len(temp_pept)) * 20)      ## Penalty for gaps
            for i, (aa, bb) in enumerate(zip(pept[:min_len], temp_pept[:min_len])):
                if i == anch_1 and aa == bb:           # Identity prize 
                    score += 3
                elif i == anch_2 and aa == bb:         # Identity prize
                    score += 1
                try:
                    gain = MatrixInfo.pam30[aa, bb]
                    if gain > 0:
                        gain += 2                      #This is summed on the top of the identity prize abobe
                    score += gain
                except KeyError:
                    try:
                        gain = MatrixInfo.pam30[bb, aa]
                        if gain > 0:
                            gain += 2                      #This is summed on the top of the identity prize abobe
                        score += gain
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
        non_modelled.append((target_id, pept, allele, "NA", "NA", "NA"))
        continue
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
    
    outdir = ('outputs/%s/%s_%s' %(outdir_name, template_ID.lower(), target_id))
    try:
        os.mkdir(outdir)
    except FileExistsError:
        print('WARNING: Existing directory.')
        continue
    os.chdir(outdir)
    
    #Preparing template and target sequences for the .ali file, launching Muscle

    template_seqr = SeqRecord(Seq(sequences[0]['M'], IUPAC.protein), id=template_ID, name = ID)
    target_seqr = SeqRecord(Seq(M_chain, IUPAC.protein), id=target_id, name = target_id)
    #target_seqr = SeqRecord(Seq(sequences[0]['M'], IUPAC.protein), id=target_id, name =target_id) #FAKE NAME
    SeqIO.write((template_seqr, target_seqr), "%s.fasta" %template_ID, "fasta")
    
    os.system('muscle -in %s.fasta -out %s.afa -quiet' %(template_ID, template_ID))
    os.system('rm %s.fasta' %template_ID)

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
            final_alifile.write('structure:../../../data/PDBs/%s_MP.pdb:%s:M:%s:P::::\n' %(template_ID, str(sequences[0]['M_st_ID']), str(len(template_pept))))
            template_ini_flag = True
            i += 1
        elif line.startswith('>') and i == 1:
            final_alifile.write('/' + template_pept + '*')
            final_alifile.write('\n')
            final_alifile.write('\n>P1;' + line.split(':')[0].strip('>'))
            final_alifile.write('sequence:::::::::\n')
            target_ini_flag = True
        else:
            final_alifile.write(line.rstrip())
    final_alifile.write('/' + str(pept) + '*')
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
                modscript.write(line %(template_ID, target_id))
            else:
                modscript.write(line)
        cmd_m_temp.close()
    
    os.popen('python3 cmd_modeller_ini.py').read()
    
    if remove_temp_outputs:
        os.system('rm cmd_modeller_ini.py')
    
    # Calculating all Atom contacts
    if "contact-chainID_allAtoms" not in os.listdir('../../../modelling_scripts'):
        os.popen('g++ ../../../modelling_scripts/contact-chainID_allAtoms.cpp -o ../../../modelling_scripts/contact-chainID_allAtoms').read()
    os.popen('../../../modelling_scripts/contact-chainID_allAtoms %s.ini %s > all_contacts_%s.list' %(target_id, cutoff, template_ID)).read()
    
    #Selecting only the anchors contacts
    
    real_anchor_2 = None
    anch_1_same = False
    anch_2_same = False
    
    if template_pept[anch_1] == pept[anch_1]:
        anch_1_same = True
    if (anch_2 + 1) == len(template_pept):
        if template_pept[anch_2] == pept[anch_2]:
            anch_2_same = True
    else:
        if template_pept[-1] == pept[anch_2]:
            anch_2_same = True
    
    if anch_2 >= len(template_pept):
        real_anchor_2 = len(template_pept)  #????
    
    ### Writing anchors contact list ###
    
    with open( 'all_contacts_%s.list' %template_ID, 'r') as contacts:                             # Template contacts
        with open('contacts_%s.list' %template_ID, 'w') as output:
            if real_anchor_2:                                                                     ### If the target peptide is longer than the templtate peptide
                for line in contacts:
                    #print(line[30:33])
                    p_aa_id = line.split("\t")[7]                                                 ### position id of the template peptide residue
                    p_atom = line.split("\t")[8]                                                  ### atom name of the template peptide residue
                    m_aa_id = (line.split("\t")[2]).split(' ')[0]
                    if anch_1_same == True:                                                       ### If the target anchor 1 residue is the same as the template anchor 1 residue
                        if int(p_aa_id) == anch_1:
                            output.write(line)
                    else:
                        if int(p_aa_id) == anch_1 and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line)
                    if anch_2_same == True:                                                       ### If the target anchor 2 residue is the same as the template anchor 2 residue
                        if int(p_aa_id) == real_anchor_2:
                            output.write(line[:30] + str(anch_2) + line[34:])
                    else:
                        if int(p_aa_id) == real_anchor_2 and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line[:30] + str(anch_2) + line[34:])
            else:
                for line in contacts:
                    #print(line.split("\t"))
                    p_aa_id = line.split("\t")[7]
                    p_atom = line.split("\t")[8]
                    m_aa_id = (line.split("\t")[2]).split(' ')[0]
                    if anch_1_same == True:                                                       ### If the target anchor 1 residue is the same as the template anchor 1 residue
                        if int(p_aa_id) == anch_1:
                            output.write(line)
                    else:
                        if int(p_aa_id) == anch_1 and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line)
                    if anch_2_same == True:                                                       ### If the target anchor 2 residue is the same as the template anchor 2 residue
                        if int(p_aa_id) == anch_2:
                            output.write(line)
                    else:
                        if int(p_aa_id) == anch_2 and ('CA' in p_atom or 'CB' in p_atom):
                            output.write(line)
                    #if (int(p_aa_id) == anch_1 or int(p_aa_id) == anch_2) and ('CA' in p_atom or 'CB' in p_atom):
                    #        output.write(line)
    
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
    
    
    with open('cmd_modeller.py', 'w') as modscript:
        cmd_m_temp = open('../../../modelling_scripts/cmd_modeller_template.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line %final_alifile_name)
            elif 'knowns' in line:
                modscript.write(line %(template_ID, target_id))
            else:
                modscript.write(line)
        cmd_m_temp.close()
    
    smt = time.time()
    os.popen('python3 cmd_modeller.py > modeller.log').read()
    
    t2 = time.time()
    tf = t2 - t1
    mt = t2 - smt
    print('The modelling took %i seconds' %mt)

    if mt < 3:
        non_modelled.append((target_id, pept, allele, template_ID, template_pept, IDD[template_ID]['allele']))
    else:
        ########################
        # BENCHMARKING #
        #TODO: 
        # [V] extract molpdf score, 
        # [0]  (optional) match all models with target .pdb, 
        # [V]  calculate lRMSD,
        # []  add molpdf and RMSD to a .tsv file
        
        #extracting molpdf and DOPE scores from modeller.log
        os.system('python ../../../modelling_scripts/get_molpdf_dope_scores.py modeller.log')
        
        #Calculatin RMSD with target real structure
        
        #os.system('python ../../../tools/make_file_lists_for_rmsd.py')
        os.system('cp ../../../data/PDBs/%s_MP.pdb ./' %target_id)
        
        os.popen('pdb_reres -1 %s_MP.pdb > reres_%s.pdb' %(target_id, target_id)).read()
        os.popen('bash ../../../tools/map_2_pdb.sh %s.BL00010001.pdb reres_%s.pdb > ref.pdb' %(target_id, target_id)).read()
        os.popen('python ../../../tools/pdb_fast_lzone_mhc.py ref.pdb').read()
        
        
        #TODO: Add here Haddock protocol on structures matching
        with open('file.list', 'w') as out:
            for f in os.listdir('./'):
                if f.startswith(target_id+'.'+filename_start) and f.endswith(filename_end):
                    out.write('matched_'+ f + '\n')
                    os.popen('bash ../../../tools/map_2_pdb.sh ref.pdb %s >  matched_%s' %(f,f)).read()
        
        os.popen('../../../tools/CA_l-rmsd-calc.csh ref file.list').read()
        os.popen('../../../tools/CA_backbone_l-rmsd-calc.csh ref file.list').read()
        os.popen('../../../tools/CA_backbone_CB_l-rmsd-calc.csh ref file.list').read()
        

        #os.system('python ../../../tools/pdb_lzone_2files.py %s_MP.pdb %s.BL00010001.pdb' %(target_id, target_id))
        #os.system('mv ref.lzone %s_MP.lzone' %target_id)
        #os.system('../../../tools/l-rmsd-calc.csh %s_MP file.list' %target_id)
        #os.system('unlink %s_MP.pdb' %target_id)
        
        models_dict = {}
        header = ['Model', 'molpdf', 'DOPE']
        with open('molpdf_DOPE.tsv', 'r') as molfile:
            spamreader = csv.reader(molfile, delimiter='\t')
            for i, row in enumerate(spamreader):
                if i != 0:
                    models_dict['matched_'+row[0]] = [float(row[1]), float(row[2])]
                    
        with open('./CA-l-RMSD.dat') as cafile:
            r = csv.reader(cafile, delimiter=' ')
            header.append('CA_l-RMSD')
            for row in r:
                models_dict[row[0]].append(float(row[1]))
        
        with open('./BB-l-RMSD.dat') as cafile:
            r = csv.reader(cafile, delimiter=' ')
            header.append('BB_l-RMSD')
            for row in r:
                models_dict[row[0]].append(float(row[1]))
        
        with open('./BB-CB-l-RMSD.dat') as cafile:
            r = csv.reader(cafile, delimiter=' ')
            header.append('BB_CB_l-RMSD')
            for row in r:
                try:
                    models_dict[row[0]].append(float(row[1]))
                except:
                    pass
        '''
        with open('l-RMSD.dat', 'r') as rmsdfile:
            for line in rmsdfile:
                row = line.split(' ')
                models_dict[row[0]].append(float(row[1]))

        
        with open('final_scores.tsv', 'wt') as scoreout:
            tsv_writer = csv.writer(scoreout, delimiter='\t')
            tsv_writer.writerow(['Model', 'molpdf', 'DOPE', 'l-RMSD'])
            for model in models_dict:
                tsv_writer.writerow([model, models_dict[model][0], models_dict[model][1], models_dict[model][2]])
        '''
        
        with open('./rmsds_and_final_scores.tsv', 'wt') as outfile:
            tw = csv.writer(outfile, delimiter='\t')
            tw.writerow(header)
            for key in models_dict:
                try:
                    tw.writerow([key, models_dict[key][0], models_dict[key][1], models_dict[key][2],
                                  models_dict[key][3], models_dict[key][4]])
                except:
                    tw.writerow([key, models_dict[key][0], models_dict[key][1], models_dict[key][2],
                                  models_dict[key][3], 'N/A'])
    
        '''
        molsort = sorted(models_dict.items(), key=lambda x:x[1][0])
        bestmol = sum(x[1][2] for x in molsort[0:5])/5.0
        dopesort = sorted(models_dict.items(), key=lambda x:x[1][1])
        bestdope = sum(x[1][2] for x in dopesort[0:5])/5.0
        best_rmsds.append([target_id, bestmol, bestdope])
        '''
        molsort = sorted(models_dict.items(), key=lambda x:x[1][0])
        bb_bestmol = sum(x[1][2] for x in molsort[0:5])/5.0
        ca_bestmol = sum(x[1][3] for x in molsort[0:5])/5.0
        try:
            ca_bb_cb_bestmol = sum(x[1][4] for x in molsort[0:5])/5.0
        except:
            ca_bb_cb_bestmol = 'N/A'
        
        dopesort = sorted(models_dict.items(), key=lambda x:x[1][1])
        bb_bestdope = sum(x[1][2] for x in dopesort[0:5])/5.0
        ca_bestdope = sum(x[1][3] for x in dopesort[0:5])/5.0
        try:
            ca_bb_cb_bestdope = sum(x[1][4] for x in dopesort[0:5])/5.0
        except:
            ca_bb_cb_bestdope = 'N/A'

        best_rmsds[target_id] = [ca_bestmol, bb_bestmol, ca_bb_cb_bestmol, ca_bestdope, bb_bestdope, ca_bb_cb_bestdope, tf]

    os.chdir('../../../')
    

dictfile = open("outputs/%s/all_best_RMSDs.pkl" %outdir_name, "wb")
pickle.dump(best_rmsds, dictfile)
dictfile.close()

with open("outputs/%s/all_best_RMSDs.tsv" %outdir_name, 'wt') as finalfile:
    tw = csv.writer(finalfile, delimiter='\t')
    tw.writerow(['Target ID', 'Top 5 molpdf CA RMDS','Top 5 molpdf BB RMDS',
             'Top 5 molpdf BB_CB RMDS', 'Top 5 DOPE CA RMDS',
             'Top 5 DOPE BB RMDS', 'Top 5 DOPE BB_CB RMDS', 'Modelling Time'])
    for key in best_rmsds:
        tw.writerow([key, best_rmsds[key][0], best_rmsds[key][1], best_rmsds[key][2],
                         best_rmsds[key][3], best_rmsds[key][4], best_rmsds[key][5], best_rmsds[key][6]])
    

with open('outputs/%s/non_modelled.csv' %outdir_name, 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['QUERY', 'NEOANTIGEN', 'QUERY ALLELE', 'TEMPLATE ID', 'TEMPLATE PEPTIDE', 'TEMPLATE ALLELE'])
    for query in non_modelled:
        tsv_writer.writerow(query)
        
'''
with open('outputs/%s/best_RMSDs.tsv' %outdir_name, 'wt') as outfile:
    tw = csv.writer(outfile, delimiter='\t')
    tw.writerow(['Target ID', 'Top 5 molpdf RMDS', 'Top 5 DOPE RMDS'])
    for query in best_rmsds:
        tw.writerow(query)
'''

final_time = time.time()
total_time = final_time - start_time
print('The whole benchmark took: ', total_time)
