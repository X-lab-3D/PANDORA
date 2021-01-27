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
from Bio.Alphabet import IUPAC
from Bio.SubsMat import MatrixInfo

from random import choice
from joblib import Parallel, delayed
from multiprocessing import Manager

from parsing import structures_parser
from parsing import url_protocols
from parsing import utils
from parsing.get_pMHC_anchors import get_anchors_pMHCI as get_anchors
from modelling import modelling

### Retriving Dictionary with PDB IDs and chain lengths ###

url_protocols.download_unzip_imgt_structures(del_inn_files = True, del_kabat_files = True)
IDs_list = url_protocols.download_ids_imgt('MH1', out_tsv='all_MH1_IDs.tsv')
IDs_list = []
with open('PANDORA_files/data/csv_pkl_files/all_MH1_IDs.tsv', 'r') as infile:
    next(infile)
    for line in infile:
        IDs_list.append(line.replace('\n',''))
IDs_dict, bad_IDs = structures_parser.parse_pMHCI_pdbs(IDs_list)


start_time = time.time()
outdir_name = sys.argv[1]

###########################################
### Organizing dicts ###

IDD_file = open('PANDORA_files/data/csv_pkl_files/IDs_and_bad_IDs_dict.pkl', 'rb')
#IDD_file = open('PANDORA_files/data/csv_pkl_files/fake_db.pkl', 'rb')
main_IDD = pickle.load(IDD_file)
bad_IDs = pickle.load(IDD_file)
IDD_file.close()

###########################################
###  Setting variables

# tsv with: Target ID, pept seq, allele
# pept_seqs = structures_parser.get_peptides_from_csv('benchmark/cross_validation_set.tsv', 3, 4, '\t') ### This can be the input requested for the benchmark function

os.system('python benchmark/get_cv_set.py')

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
            pept_seqs.append((target_id, seq, allele))

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

remove_temp_outputs = False
cutoff = 5

maxs = []
maxsl = []
non_modelled = []
print_results = False

filename_start = 'BL00'
filename_end = '0001.pdb'

manager = Manager()
best_rmsds = manager.dict()

###########################################

def na_model(k, pept_seq, best_rmsds):
    t1 = time.time()
    #for k, pept_seq in zip(range(pepts_start, pepts_end), pept_seqs[pepts_start:pepts_end]):
    #for k, pept_seq in enumerate(pept_seqs[575:576]):

    query = 'query_' + str(k + 1)
    target_id = pept_seq[0]
    pept = pept_seq[1]
    allele = pept_seq[2]
    length = len(pept)
    try:
        M_chain = utils.get_seqs('./PANDORA_files/data/PDBs/pMHCI/%s_MP.pdb' %target_id)['M']
    except FileNotFoundError:
        print('Error in retrieving M_chain, could not find PDB')
        print('CWD: %s' %os.getcwd())
        return (target_id, pept, allele, "NA", "NA", "NA", "Invalid Peptide Length")

    if length < 8 or length > 12:
        print('###')
        print('Invalid Peptide Length. Exiting here.')
        print('###')
        return (target_id, pept, allele, "NA", "NA", "NA", "Invalid Peptide Length")

    IDD = copy.deepcopy(main_IDD)
    del IDD[target_id]

    ### Organizing Allele IDs in a dictionary ###
    allele_list = []
    for key in IDD:
        allele_list += IDD[key]['allele']
        allele_list = list(set(allele_list))

    allele_ID = {i : [] for i in allele_list}
    for key in IDD:
        for multi_allele in IDD[key]['allele']:
            allele_ID[multi_allele].append(key)

    #Check if this query has been modelled already
    ext_flag = False
    print( '## Wokring on %s, structure %s ##' %(query, target_id))
    ext_flag = ext_flag = modelling.check_model_existance(k, outdir_name, filename_start, filename_end)
    if ext_flag == True:
        print('Exiting here. This is only a DEBUG message')
        return None

    #Correct allele name depending on available alleles
    allele = modelling.allele_name_adapter(allele, allele_ID)

    #Select template
    sel_template = modelling.select_template(IDD, allele, homolog_allele, length, pept, print_results)

    ###################################
    #   CHANGING WORKING DIRECTORY
    ###################################

    sequences, empty_seqs = utils.get_pdb_seq([template_ID])

    outdir = ('outputs/%s/%s_%s' %(outdir_name, template_ID.lower(), target_id))
    try:
        os.mkdir(outdir)
    except FileExistsError:
        print('WARNING: Existing directory.')
        return None
    os.chdir(outdir)

    #Obtain template anchor positions
    try:
        anch_1, anch_2 = get_anchors('../../../data/PDBs/pMHCI/%s_MP.pdb' %target_id, rm_outfile = True)
        anch_1 -= 1
        anch_2 -= 1
    except:
        os.chdir('../../../')
        return (target_id, pept, allele, template_ID, template_pept, IDD[template_ID]['allele'], 'Something gone wrong in defining target anchors')

    try:
        temp_anch_1, temp_anch_2 = get_anchors('../../../data/PDBs/pMHCI/%s_MP.pdb' %template_ID, rm_outfile = True)
        temp_anch_1 -= 1
        temp_anch_2 -= 1
    except:
        os.chdir('../../../')
        return (target_id, pept, allele, template_ID, template_pept, IDD[template_ID]['allele'], 'Something gone wrong in defining target anchors')

    #Write alignment file.
    final_alifile_name = modelling.make_ali_files(sequences[0], M_chain, template_ID, query, template_pept, pept, length, remove_temp_outputs)

    #Write and running MODELLER scripts to produce .ini file
    modelling.write_ini_scripts(anch_1, anch_2, template_ID, final_alifile_name, query)
    #If default python is not python3, change the following line
    os.popen('python cmd_modeller_ini.py').read()

    if remove_temp_outputs:
        os.system('rm cmd_modeller_ini.py')

    # Calculate all Atom contacts
    os.popen('../../../../PANDORA/tools/contact-chainID_allAtoms %s.ini %s > all_contacts_%s.list' %(target_id, cutoff, template_ID)).read()

    #Select only the anchors contacts
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
    ### Writing anchors contact list ###

    modelling.write_contact_list(anch_1, anch_2, anch_1_same, anch_2_same)
    '''
    with open( 'all_contacts_%s.list' %template_ID, 'r') as contacts:                             # Template contacts
        with open('contacts_%s.list' %template_ID, 'w') as output:                                                                ### If the target peptide is longer than the templtate peptide #TODO: This can be removed?
            for line in contacts:
                p_aa_id = line.split("\t")[7]                                                 ### position id of the template peptide residue
                p_atom = line.split("\t")[8]                                                  ### atom name of the template peptide residue
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
    '''
    if remove_temp_outputs:
        os.system('rm all_contacts_%s.list' %template_ID)

    with open('MyLoop.py', 'w') as myloopscript:
        MyL_temp = open('../../../../PANDORA/modelling/MyLoop_template.py', 'r')
        for line in MyL_temp:
            if 'self.residue_range' in line:
                myloopscript.write(line %(anch_1 + 2, anch_2))
            elif 'contact_file = open' in line:
                myloopscript.write(line %template_ID)
            else:
                myloopscript.write(line)
        MyL_temp.close()

    with open('cmd_modeller.py', 'w') as modscript:
        cmd_m_temp = open('../../../../PANDORA/modelling/cmd_modeller_template.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line %final_alifile_name)
            elif 'knowns' in line:
                modscript.write(line %(template_ID, target_id))
            else:
                modscript.write(line)
        cmd_m_temp.close()

    #Running MODELLER
    smt = time.time()
    os.popen('python3 cmd_modeller.py > modeller.log').read()

    t2 = time.time()
    tf = t2 - t1
    mt = t2 - smt
    print('The modelling took %i seconds' %mt)

    if mt < 3:
        os.chdir('../../../')
        return (target_id, pept, allele, template_ID, template_pept, IDD[template_ID]['allele'], 'Something gone wrong in the modelling')
    else:
        ########################
        # BENCHMARKING #

        #extracting molpdf and DOPE scores from modeller.log
        os.system('python ../../../parsing/get_molpdf_dope_scores.py modeller.log')

        #Calculating RMSD with target real structure

        #os.system('python ../../../tools/make_file_lists_for_rmsd.py')
        os.system('cp ../../../data/PDBs/pMHCI/%s_MP.pdb ./' %target_id)

        os.popen('pdb_reres -1 %s_MP.pdb > reres_%s.pdb' %(target_id, target_id)).read()

        for f in os.listdir('./'):
            if f.startswith(target_id+'.'+filename_start) and f.endswith(filename_end):
                break
        os.popen('bash ../../../tools/map_2_pdb.sh %s reres_%s.pdb > ref.pdb' %(f, target_id)).read()
        os.popen('python ../../../tools/pdb_fast_lzone_mhc_fitG.py ref.pdb').read()

        with open('file.list', 'w') as out:
            for f in os.listdir('./'):
                if f.startswith(target_id+'.'+filename_start) and f.endswith(filename_end):
                    out.write('matched_'+ f + '\n')
                    os.popen('bash ../../../tools/map_2_pdb.sh ref.pdb %s >  matched_%s' %(f,f)).read()

        os.popen('../../../tools/CA_l-rmsd-calc.csh ref file.list').read()
        os.popen('../../../tools/CA_backbone_l-rmsd-calc.csh ref file.list').read()
        os.popen('../../../tools/CA_backbone_CB_l-rmsd-calc.csh ref file.list').read()

        try:
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

        except:
            print('Warning: something when wrong in parsing RMSD ouputs of target %s . Moving on.' %target_id)

    os.chdir('../../../')
    return None

######################################################################################################
######################################################################################################

non_modelled = Parallel(n_jobs = n_cores)(delayed(na_model[:100])(k, pept_seq, best_rmsds) for k, pept_seq in enumerate(pept_seqs[pepts_start:pepts_end]))
non_modelled = filter(None, non_modelled)

best_rmsds = dict(best_rmsds)
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
    tsv_writer.writerow(['QUERY', 'NEOANTIGEN', 'QUERY ALLELE', 'TEMPLATE ID', 'TEMPLATE PEPTIDE', 'TEMPLATE ALLELE', 'ERROR'])
    for query in non_modelled:
        tsv_writer.writerow(query)

final_time = time.time()
total_time = final_time - start_time
print('The whole benchmark took: ', total_time)
