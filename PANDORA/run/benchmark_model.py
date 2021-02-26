#!/usr/bin/python

###    ###

import os
import sys
import time
import copy
import pickle
import csv
#import subprocess

from joblib import Parallel, delayed
from multiprocessing import Manager

sys.path.append('./')

import PANDORA
from PANDORA.Database import utils
from PANDORA.junk.parsing import get_anchors_pMHCI as get_anchors
from PANDORA.junk.modelling import modelling

### Retriving Dictionary with PDB IDs and chain lengths ###

#url_protocols.download_unzip_imgt_structures(del_inn_files = True, del_kabat_files = True)
#IDs_list = url_protocols.download_ids_imgt('MH1', out_tsv='all_MH1_IDs.tsv')

'''
IDs_list = []
with open( PANDORA.PANDORA_data + '/csv_pkl_files/all_MH1_IDs.tsv', 'r') as infile:
    next(infile)
    for line in infile:
        IDs_list.append(line.replace('\n',''))
'''

#main_IDD, bad_IDs = structures_parser.parse_pMHCI_pdbs(IDs_list)


start_time = time.time()
outdir_name = sys.argv[1]

###########################################
### Organizing dicts ###

IDD_file = open( PANDORA.PANDORA_data + '/csv_pkl_files/IDs_and_bad_IDs_dict.pkl', 'rb')
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

### FOR SHORT TESTS
#pept_seqs = pept_seqs[:40]

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
cutoff = str(5)

maxs = []
maxsl = []
non_modelled = []
print_results = True

filename_start = 'BL00'
filename_end = '0001.pdb'

manager = Manager()
best_rmsds = manager.dict()

###########################################

def na_model(k, target_info, best_rmsds):
    '''
    Prepares and performs the model for one case with the following steps:
    1. Data parsing and Preparation
    2. Template selection
    3. Anchor selection
    4. Write alingment files and MODELLER .ini file
    5. Retrive distance restrains
    6. Prepare MODELLER modelling scripts
    7. Run modeller
    8. Models benchmark (only for cross-validation purposes)

    Args:
        k(int) : query index for batch jobs. It is used to keep track of the query
        target_info(list) : list containing target PDB ID, peptide sequence and MHC allele
        best_rmsds(dict) : shared dictionary of top models l-RMSDs
    '''

    ###############################
    ########## Step 1 ############
    ###############################

    t1 = time.time()

    ### Set some useful variables
    query = 'query_' + str(k + 1)
    target_id = target_info[0]
    pept = target_info[1]
    allele = target_info[2]
    length = len(pept)

    ### Retrieve the target M chain sequence
    try:
        M_chain = utils.get_seqs(PANDORA.PANDORA_data + '/PDBs/pMHCI/%s_MP.pdb' % target_id)['M']
    except FileNotFoundError:
        print('Error in retrieving M_chain, could not find PDB')
        print('CWD: %s' %os.getcwd())
        return (target_id, pept, allele, "NA", "NA", "NA", "Invalid Peptide Length")

    if length < 8 or length > 15:
        print('###')
        print('Invalid Peptide Length. Exiting here.')
        print('###')
        return (target_id, pept, allele, "NA", "NA", "NA", "Invalid Peptide Length")


    ### Copying the main ID Dictionary to a fake one for cross validation.
    ### This IDD will have all the IDs and structure info except the target one.
    IDD = copy.deepcopy(main_IDD)
    del IDD[target_id]


    ### Organizing Allele IDs in a dictionary ###
    ### This dictionary contains alleles as keys and structure IDs as values,
    ### reporting which structures belong to which allele.
    ### Please note that since some structure can be identified as belonging to
    ### multiple alleles, thi disctionary will be ambiguous
    allele_list = []
    for key in IDD:
        allele_list += IDD[key]['allele']
        allele_list = list(set(allele_list))

    allele_ID = {i : [] for i in allele_list}
    for key in IDD:
        for multi_allele in IDD[key]['allele']:
            allele_ID[multi_allele].append(key)

    ### Check if this query has been modelled already

    ext_flag = False
    print( '## Wokring on %s, structure %s ##' %(query, target_id))
    ext_flag = modelling.check_model_existance(k, outdir_name, filename_start, filename_end)
    if ext_flag == True:
        print('Exiting here. This is only a DEBUG message')
        return None

    ### Correct allele name depending on available template alleles
    allele = modelling.allele_name_adapter(allele, allele_ID)

    ###############################
    ########## Step 2 ############
    ###############################

    ### Select template
    sel_template = modelling.select_template(IDD, target_id, allele, length, pept, print_results) # ,homolog_allele
    template_pept = sel_template[1]
    template_ID = sel_template[2]


    ### Retrieve template PDB sequences
    try:
        sequences, empty_seqs = utils.get_pdb_seq([template_ID])
    except FileNotFoundError:
        return (target_id, pept, allele, template_ID, template_pept, 'NA', 'Something gone wrong in retrieving template sequence')


    ### Set output directory
    outdir = ('PANDORA_files/outputs/%s/%s_%s' %(outdir_name, template_ID.lower(), target_id))
    try:
        os.mkdir(outdir)
    except FileNotFoundError: #FileExistsError:
        os.mkdir('./TESTDIR')
        raise Exception('TESTDIR generated')
        #print('WARNING: Existing directory.')
        #return None

    ### Changing working directory
    os.chdir(outdir)


    ###############################
    ########## Step 3 ############
    ###############################

    ### Obtain template anchor positions
    try:
        anch_1, anch_2 = get_anchors( PANDORA.PANDORA_data + '/PDBs/pMHCI/%s_MP.pdb' %target_id, rm_outfile = False)
        anch_1 -= 1
        anch_2 -= 1
    except IndexError:
        print(target_id)
        os.chdir('../../../../')
        return (target_id, pept, allele, template_ID, template_pept, IDD[template_ID]['allele'], 'Something gone wrong in defining target anchors')

    try:
        temp_anch_1, temp_anch_2 = get_anchors(PANDORA.PANDORA_data + '/PDBs/pMHCI/%s_MP.pdb' %template_ID, rm_outfile = False)
        temp_anch_1 -= 1
        temp_anch_2 -= 1
    except IndexError:
        print(template_ID)
        os.chdir('../../../../')
        return (target_id, pept, allele, template_ID, template_pept, IDD[template_ID]['allele'], 'Something gone wrong in defining template anchors')

    ###############################
    ########## Step 4 ############
    ###############################

    ### Write alignment file.
    final_alifile_name = modelling.make_ali_files(sequences[0], M_chain, template_ID, target_id,
                                                  template_pept, temp_anch_1, temp_anch_2, pept,
                                                  anch_1, anch_2, length, remove_temp_outputs)

    ### Write and running MODELLER scripts to produce .ini file
    modelling.write_ini_scripts(anch_1, anch_2, template_ID, final_alifile_name, target_id)

    ### If default python is not python3, change the following line
    os.popen('python cmd_modeller_ini.py').read()

    if remove_temp_outputs:
        os.system('rm cmd_modeller_ini.py')

    ###############################
    ########## Step 5 ############
    ###############################

    ###  Calculate all Atom contacts
    #s.popen('../../../../PANDORA/tools/contact-chainID_allAtoms %s.ini %s > all_contacts_%s.list' %(target_id, cutoff, template_ID)).read()

    dist_file ='all_contacts_'+ template_ID +'.list'
    os.system(PANDORA.PANDORA_path + '/tools/contact-chainID_allAtoms '
              + target_id + '.ini ' + cutoff +' > ' + dist_file)

    ### Writing anchors contact list ###

    modelling.write_contact_list(pept, anch_1, anch_2, template_pept, temp_anch_1, temp_anch_2, template_ID)

    if remove_temp_outputs:
        os.system('rm all_contacts_%s.list' %template_ID)

    ###############################
    ########## Step 6 ############
    ###############################

    modelling.write_modelling_scripts(anch_1, anch_2, template_ID, target_id, final_alifile_name)

    ###############################
    ########## Step 7 ############
    ###############################

    #Running MODELLER
    smt = time.time()
    os.popen('python3 cmd_modeller.py > modeller.log').read()

    t2 = time.time()
    tf = t2 - t1
    mt = t2 - smt
    print('The modelling took %i seconds' %mt)

    if mt < 3:
        os.chdir('../../../../')
        return (target_id, pept, allele, template_ID, template_pept, IDD[template_ID]['allele'], 'Something gone wrong in the modelling')
    else:
        ###############################
        ########## Step 1 ############
        ###############################

        ########################
        # BENCHMARKING #
        ########################

        #extracting molpdf and DOPE scores from modeller.log
        os.system('python ' + PANDORA.PANDORA_path + '/parsing/get_molpdf_dope_scores.py modeller.log')

        #Calculating RMSD with target real structure

        #os.system('python ../../../tools/make_file_lists_for_rmsd.py')
        os.system('cp ' + PANDORA.PANDORA_data +'/PDBs/pMHCI/%s_MP.pdb ./' %target_id)

        os.popen('pdb_reres -1 %s_MP.pdb > reres_%s.pdb' %(target_id, target_id)).read()

        for f in os.listdir('./'):
            if f.startswith(target_id+'.'+filename_start) and f.endswith(filename_end):
                break
        #os.popen('bash ' + PANDORA.PANDORA_path + '/tools/map_2_pdb.sh %s reres_%s.pdb > ref.pdb' %(f, target_id)).read()
        os.popen('python ' + PANDORA.PANDORA_path + '/tools/map_2_pdb.py %s reres_%s.pdb > ref.pdb' %(f, target_id)).read()
        os.popen('python ' + PANDORA.PANDORA_path + '/tools/pdb_fast_lzone_mhc_fitG.py ref.pdb').read()

        with open('file.list', 'w') as out:
            for f in os.listdir('./'):
                if f.startswith(target_id+'.'+filename_start) and f.endswith(filename_end):
                    out.write('matched_'+ f + '\n')
                    os.popen('bash ' + PANDORA.PANDORA_path + '/tools/map_2_pdb.sh ref.pdb %s >  matched_%s' %(f,f)).read()

        os.popen( PANDORA.PANDORA_path + '/tools/CA_l-rmsd-calc.csh ref file.list').read()
        os.popen(PANDORA.PANDORA_path + '/tools/CA_backbone_l-rmsd-calc.csh ref file.list').read()
        os.popen(PANDORA.PANDORA_path + '/tools/CA_backbone_CB_l-rmsd-calc.csh ref file.list').read()

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

    os.chdir('../../../../')
    return None

######################################################################################################
######################################################################################################

non_modelled = Parallel(n_jobs = n_cores, verbose = 1)(delayed(na_model)(k, target_info, best_rmsds) for k, target_info in enumerate(pept_seqs[pepts_start:pepts_end]))
non_modelled = filter(None, non_modelled)

best_rmsds = dict(best_rmsds)
dictfile = open("PANDORA_files/outputs/%s/all_best_RMSDs.pkl" %outdir_name, "wb")
pickle.dump(best_rmsds, dictfile)
dictfile.close()

with open("PANDORA_files/outputs/%s/all_best_RMSDs.tsv" %outdir_name, 'wt') as finalfile:
    tw = csv.writer(finalfile, delimiter='\t')
    tw.writerow(['Target ID', 'Top 5 molpdf CA RMDS','Top 5 molpdf BB RMDS',
             'Top 5 molpdf BB_CB RMDS', 'Top 5 DOPE CA RMDS',
             'Top 5 DOPE BB RMDS', 'Top 5 DOPE BB_CB RMDS', 'Modelling Time'])
    for key in best_rmsds:
        tw.writerow([key, best_rmsds[key][0], best_rmsds[key][1], best_rmsds[key][2],
                         best_rmsds[key][3], best_rmsds[key][4], best_rmsds[key][5], best_rmsds[key][6]])


with open('PANDORA_files/outputs/%s/non_modelled.csv' %outdir_name, 'wt') as outfile:
    tsv_writer = csv.writer(outfile, delimiter='\t')
    tsv_writer.writerow(['QUERY', 'NEOANTIGEN', 'QUERY ALLELE', 'TEMPLATE ID', 'TEMPLATE PEPTIDE', 'TEMPLATE ALLELE', 'ERROR'])
    for query in non_modelled:
        tsv_writer.writerow(query)

final_time = time.time()
total_time = final_time - start_time
print('The whole benchmark took: ', total_time)
