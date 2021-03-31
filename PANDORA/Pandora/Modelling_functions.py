from Bio.Align import substitution_matrices
import os
import PANDORA
import dill
import pickle
from PANDORA.PMHC import Model
# from Bio import Align
from PANDORA.Pandora import Align

def check_target_template(target, template):
    """ Checks if the target and the template are the same. If the user gave sequence info in the target, use that, else
        use the allele type.

    Args:
        target: (:obj:`Target`): Target object
        template: (:obj:`Template`): Template object

    Returns: (bool): True if target/template are the same, False if they are not.

    """
    out = False
    # Check if target peptide and template peptide are the same
    if target.peptide == template.peptide:
        # If the target has no sequence information, use allele type
        if target.M_chain_seq == '':
            # Check if the allele of target and template are the same
            if any(x in template.allele_type for x in target.allele_type):
                out = True

        # If the target has sequence information..
        elif target.M_chain_seq != '':
            # For MHCI, check if the M chain sequence of target and template are the same
            if target.MHC_class == 'I':
                if target.M_chain_seq == template.M_chain_seq:
                    out = True
            # For MHCII, check if the M and N chain sequence of target and template are the same
            elif target.MHC_class == 'II' and target.N_chain_seq != '':
                if target.M_chain_seq == template.M_chain_seq and target.N_chain_seq == template.N_chain_seq:
                    out = True
    if out:
        print('\n\t---- THE TARGET HAS THE SAME PEPTIDE AND ALLELE/SEQUENCE INFORMATION AS THE TEMPLATE ----')
        print('\tYou can find it at: http://www.imgt.org/3Dstructure-DB/cgi/details.cgi?pdbcode=%s\n' %(template.id))

    return out






def check_presence(target, database, seq_based_templ_selection = False):
    ''' Checks if the target the user submitted, already exists in has a template in the database with the same allele
        and peptide.

    Args:
        target: Target object
        database: Database object
        seq_based_templ_selection: bool, select the template based on the chain sequences.

    Returns: bool/Template object. If the target is already in the db, return the Template, otherwise return False
    
    '''
    putative_templates = []
    target_in_db = False
    if not seq_based_templ_selection:
        # For MHC I
        if target.MHC_class == 'I':
            # Check if there are templates with the same alleles
            for id in database.MHCI_data:
                if any(x in database.MHCI_data[id].allele_type for x in target.allele_type):
                    putative_templates.append(id)
            # Check if there is a putative template that also has the same peptide as the target
            for i in putative_templates:
                if database.MHCI_data[i].peptide == target.peptide:
                    target_in_db = database.MHCI_data[i]
        # For MHC II
        elif target.MHC_class == 'II':
            # Check if there are templates with the same alleles
            for id in database.MHCII_data:
                if any(x in database.MHCII_data[id].allele_type for x in target.allele_type):
                    putative_templates.append(id)
            # Check if there is a putative template that also has the same peptide as the target
            for i in putative_templates:
                if database.MHCII_data[i].peptide == target.peptide:
                    target_in_db = database.MHCII_data[i]

    elif seq_based_templ_selection:
        # Check for MHC I
        if target.MHC_class == 'I':
            # Check if there are templates with the same M chain sequence
            for id in database.MHCI_data:
                if database.MHCI_data[id].M_chain_seq == target.M_chain_seq:
                    putative_templates.append(id)
            # Check if there is a putative template that also has the same peptide as the target
            for i in putative_templates:
                if database.MHCI_data[i].peptide == target.peptide:
                    target_in_db = database.MHCI_data[i]
        # Check for MHC I
        if target.MHC_class == 'II':
            # Check if there are templates with the same M chain sequence
            for id in database.MHCII_data:
                if database.MHCII_data[id].M_chain_seq == target.M_chain_seq:
                    if database.MHCII_data[id].N_chain_seq == target.N_chain_seq:
                        putative_templates.append(id)
            # Check if there is a putative template that also has the same peptide as the target
            for i in putative_templates:
                if database.MHCII_data[i].peptide == target.peptide:
                    target_in_db = database.MHCII_data[i]

    return target_in_db


def find_template(target, database, seq_based_templ_selection = False, benchmark=False):
    ''' Selects the template structure that is best suited as template for homology modelling of the target

    Args:
        target: Target object
        database: Database object
        seq_based_templ_selection: (bool) Use template selection based on template sequences instead of allele.

    Returns: Template object

    '''
    # Check if the target is already in the database
    # templ_present = check_presence(target, database, seq_based_templ_selection = seq_based_templ_selection)
    # If the template is already present in the db and you're not benchmarking, return this template
    # if templ_present and not benchmark:
    #     print('\n\t---- PANDORA FOUND A TEMPLATE WITH THE SAME ALLELE AND PEPTIDE SEQUENCE AS THE TARGET ----\n')
    #     return templ_present, True

    # def find_best_template_peptide(target, templates):
    #     ''' Finds the best template based on pairwise alignment
    #
    #     Args:
    #         target: (:obj:`Target`): Target object
    #         templates: (:obj:`list` of :obj:`Template`): A list of putative templates to align and score
    #
    #     Returns: (str): ID of the template with the best matching peptide.
    #
    #     '''
    #     # Find the alignment scores of the peptides of all putative templates with the target peptide
    #     ali = Align.Align2(target=target, template=templates)
    #     ali.align_peptides()
    #
    #     # Find the peptide with the highest alignment score. If there are multiple, take the first one with same
    #     # same anchor positions
    #     template_id = max((v, k) for k, v in ali.pept_alignment_scores.items())[1]
    #     for k, v in ali.pept_alignment_scores.items():
    #         if v[0] == max((v, k) for k, v in ali.pept_alignment_scores.items())[0][0]:
    #             if v[2] == target.anchors:
    #                 template_id = k
    #
    #     return template_id


    if not seq_based_templ_selection:

        PAM30 = substitution_matrices.load('PAM30')

        ## For MHC I
        if target.MHC_class == 'I':

            # Find template structures with matching alleles
            putative_templates = {}
            for ID in database.MHCI_data:
                if benchmark:
                    if ID != target.id:
                        if any(x in database.MHCI_data[ID].allele_type for x in target.allele_type):
                            putative_templates[ID] = list(
                                set(target.allele_type) & set(database.MHCI_data[ID].allele_type))
                else:
                    if any(x in database.MHCI_data[ID].allele_type for x in target.allele_type):
                        putative_templates[ID] = list(
                            set(target.allele_type) & set(database.MHCI_data[ID].allele_type))
            # update dict with ID:all matching alleles

            # Find the putative template with the best matching peptide
            # template_id = find_best_template_peptide(target=target,
            #                                          templates=[database.MHCII_data[i] for i in putative_templates])

            pos_list = []
            for ID in putative_templates:
                score = 0
                temp_pept = database.MHCI_data[ID].peptide
                min_len = min([len(target.peptide), len(temp_pept)])
                score -= ((abs(len(target.peptide) - len(temp_pept)) ** 2.4))  # !!!  ## Gap Penalty
                for i, (aa, bb) in enumerate(zip(target.peptide[:min_len], temp_pept[:min_len])):
                    try:
                        # gain = MatrixInfo.pam30[aa, bb]
                        gain = PAM30[aa, bb]
                        score += gain
                    except KeyError:
                        try:
                            # gain = MatrixInfo.pam30[bb, aa]
                            gain = PAM30[bb, aa]
                            score += gain
                        except KeyError:
                            score = -50
                            pass
                pos_list.append((score, temp_pept, ID))

            # Take the putative template with the max scoring peptide
            template_id = pos_list[[i[0] for i in pos_list].index(max([i[0] for i in pos_list]))][2]
            # Return the Template object of the selected template that will be used for homology modelling



            return database.MHCI_data[template_id], check_target_template(target, database.MHCI_data[template_id])


        ## For MHC II
        if target.MHC_class == 'II':

            # Find template structures with matching alleles
            putative_templates = {}
            for ID in database.MHCII_data:
                if benchmark:
                    if ID != target.id:
                        if any(x in database.MHCII_data[ID].allele_type for x in target.allele_type):
                            putative_templates[ID] = list(
                                set(target.allele_type) & set(database.MHCII_data[ID].allele_type))
                    else:
                        if any(x in database.MHCII_data[ID].allele_type for x in target.allele_type):
                            putative_templates[ID] = list(
                                set(target.allele_type) & set(database.MHCII_data[ID].allele_type))

            # Find the peptide with the highest alignment score. If there are multiple, take the first one with same
            # same anchor positions

            # template_id = find_best_template_peptide(target=target,
            #                                          templates=[database.MHCII_data[i] for i in putative_templates])

            # Find the putative template with the best matching peptide
            pos_list = []
            for ID in putative_templates:
                score = 0
                temp_pept = database.MHCII_data[ID].peptide
                min_len = min([len(target.peptide), len(temp_pept)])
                score -= ((abs(len(target.peptide) - len(temp_pept)) ** 2.4))  # !!!  ## Gap Penalty
                for i, (aa, bb) in enumerate(zip(target.peptide[:min_len], temp_pept[:min_len])):
                    try:
                        # gain = MatrixInfo.pam30[aa, bb]
                        gain = PAM30[aa, bb]
                        score += gain
                    except KeyError:
                        try:
                            # gain = MatrixInfo.pam30[bb, aa]
                            gain = PAM30[bb, aa]
                            score += gain
                        except KeyError:
                            score = -50
                            pass
                pos_list.append((score, temp_pept, ID))

            # Take the putative template with the max scoring peptide
            template_id = pos_list[[i[0] for i in pos_list].index(max([i[0] for i in pos_list]))][2]
            # Return the Template object of the selected template that will be used for homology modelling

            return database.MHCII_data[template_id], check_target_template(target, database.MHCI_data[template_id])

    # Sequence based template search if the sequences of the target are provided
    elif target.M_chain_seq != '' and seq_based_templ_selection:

        if target.MHC_class == 'I':

            # define target sequences
            tar_seq = database.MHCI_data[target.id].M_chain_seq
            tar_pept = database.MHCI_data[target.id].peptide
            # keep track of alignment scores
            scores = {}
            # Perform a pairwise alignment of the target and all templates for the MHC M chain and peptide
            for i in database.MHCI_data:
                aligner = Align.PairwiseAligner()
                aligner.substitution_matrix = substitution_matrices.load("BLOSUM80")  # PAM30 for pept??

                M_score = aligner.align(tar_seq, database.MHCI_data[i].M_chain_seq).score
                P_score = aligner.align(tar_pept, database.MHCI_data[i].peptide).score

                scores[i] = (M_score, P_score)
            # Remove the target structure from this dict, you cannot select the target as template
            scores.pop(target.id, None)
            # take the 10 best scoring templates
            best_MHCs = sorted(scores, key=scores.get, reverse=True)[:10]
            # take the template with the best scoring peptide
            best_template = max((v[1], k) for k, v in scores.items() if k in best_MHCs)[1]

            return database.MHCI_data[best_template], check_target_template(target, database.MHCI_data[template_id])

        if target.MHC_class == 'II':
            # define target sequences
            tar_seq = database.MHCII_data[target.id].M_chain_seq + database.MHCII_data[target.id].N_chain_seq
            tar_pept = database.MHCII_data[target.id].peptide
            # keep track of alignment scores
            scores = {}

            for i in database.MHCII_data:
                aligner = Align.PairwiseAligner()
                aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")  # or PAM30 ??

                temp_seq = database.MHCII_data[i].M_chain_seq + database.MHCII_data[i].N_chain_seq
                MN_score = aligner.align(tar_seq, temp_seq).score
                P_score = aligner.align(tar_pept, database.MHCII_data[i].peptide).score

                scores[i] = (MN_score, P_score)
            # Remove the target structure from this dict, you cannot select the target as template
            scores.pop(target.id, None)
            # take the 10 best scoring templates
            best_MHCs = sorted(scores, key=scores.get, reverse=True)[:10]
            # take the template with the best scoring peptide
            best_template = max((v[1], k) for k, v in scores.items() if k in best_MHCs)[1]

            return database.MHCII_data[best_template], check_target_template(target, database.MHCI_data[template_id])


def write_ini_script(target, template, alignment_file, output_dir):
    ''' Writes the MyLoop.py and cmd_modeller_ini.py files. This function takes two template python scripts and fills
        in the required information: Anchor positions for the MyLoop file and structure name + alignment file for the
        cmd_modeller_ini file.

    Args:
        target: Target object
        template: Template object
        alignment_file: (string) path to alignment file
        output_dir: (string) path to output directory

    '''

    anch = target.anchors

    if target.MHC_class == 'I':
        with open(output_dir+ '/MyLoop.py', 'w') as myloopscript:
            MyL_temp = open(PANDORA.PANDORA_path + '/Pandora/MyLoop_template.py', 'r')
            for line in MyL_temp:
                if 'self.residue_range' in line:
                    myloopscript.write(line % (anch[0], anch[-1]))
                elif 'SPECIAL_RESTRAINTS_BREAK' in line:
                    break
                elif 'contact_file = open' in line:
                    myloopscript.write(line % template_ID)
                else:
                    myloopscript.write(line)
            MyL_temp.close()

    if target.MHC_class == 'II':
        with open(output_dir + '/MyLoop.py', 'w') as myloopscript:
            MyL_temp = open(PANDORA.PANDORA_path + '/Pandora/MyLoop_template_II.py', 'r')
            for line in MyL_temp:
                if 'self.residue_range' in line:
                    myloopscript.write(line % (1, anch[0])) # write the first anchor
                    for i in range(len(anch)-1): # Write all the inbetween acnhors if they are there
                        myloopscript.write(line % (anch[i] + 2, anch[i+1]))
                    myloopscript.write(line % (anch[-1] + 2, len(target.peptide))) # Write the last anchor
                elif 'SPECIAL_RESTRAINTS_BREAK' in line:
                    break
                elif 'contact_file = open' in line:
                    myloopscript.write(line % template_ID)
                else:
                    myloopscript.write(line)
            MyL_temp.close()

    with open(output_dir.replace('\\ ', ' ') + '/cmd_modeller_ini.py', 'w') as modscript:
        cmd_m_temp = open(PANDORA.PANDORA_path + '/Pandora/cmd_modeller_ini.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line % os.path.basename(alignment_file))
            elif 'knowns' in line:
                modscript.write(line % (template.id, target.id))
            else:
                modscript.write(line)
        cmd_m_temp.close()


def write_modeller_script(target, template, alignment_file, output_dir, n_models=20, stdev=0.1, n_jobs = None):
    ''' Write script that refines the loops of the peptide

    Args:
        target: Target object
        template: Template object
        alignment_file: (string) path to alignment file
        output_dir: (string) path to output directory
        n_models:  (int) number of models modeller generates per run
        stdev: (float) standard deviation of modelling restraints. Higher = more flexible restraints.

    '''



    anch = target.anchors

    if target.MHC_class == 'I':
        with open(output_dir.replace('\\ ', ' ') + '/MyLoop.py', 'w') as myloopscript:
            MyL_temp = open(PANDORA.PANDORA_path + '/Pandora/MyLoop_template.py', 'r')
            for line in MyL_temp:
                if 'self.residue_range' in line:
                    myloopscript.write(line %(anch[0], anch[-1]))  # write the first anchor
                elif 'contact_file = open' in line:
                    myloopscript.write(line %(target.id))
                elif 'STDEV MARKER' in line:
                    myloopscript.write(line %(stdev))
                else:
                    myloopscript.write(line)
            MyL_temp.close()

    if target.MHC_class == 'II':
        with open(output_dir.replace('\\ ', ' ') + '/MyLoop.py', 'w') as myloopscript:
            MyL_temp = open(PANDORA.PANDORA_path + '/Pandora/MyLoop_template_II.py', 'r')
            for line in MyL_temp:
                if 'self.residue_range' in line:
                    myloopscript.write(line % (1, anch[0]))  # write the first anchor
                    for i in range(len(anch) - 1):  # Write all the inbetween acnhors if they are there
                        myloopscript.write(line %(anch[i] + 2, anch[i + 1]))
                    myloopscript.write(line %(anch[-1] + 2, len(target.peptide)))  # Write the last anchor
                elif 'contact_file = open' in line:
                    myloopscript.write(line %(target.id))
                elif 'STDEV MARKER' in line:
                    myloopscript.write(line %(stdev))
                else:
                    myloopscript.write(line)
            MyL_temp.close()

    with open(output_dir.replace('\\ ', ' ') + '/cmd_modeller.py', 'w') as modscript:
        cmd_m_temp = open(PANDORA.PANDORA_path + '/Pandora/cmd_modeller_template.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line %(os.path.basename(alignment_file)))
            elif 'knowns' in line:
                modscript.write(line %(template.id, target.id))
            elif 'a.loop.ending_model' in line:
                modscript.write(line % (n_models))
            else:
                modscript.write(line)
        cmd_m_temp.close()


def run_modeller(output_dir, target, python_script = 'cmd_modeller.py', benchmark = False, pickle_out = True, keep_IL = False):
    ''' Perform the homology modelling.

    Args:
        output_dir: (string) path to output directory
        target: Target object
        python_script:  (string) path to script that performs the modeller modelling. cmd_modeller.py
        benchmark: (bool) Perform L-RMSD calculations? only works if the target id is an existing pdb id
        pickle_out: (bool) Save a .pkl with the results

    Returns: (list) of Model objects

    '''
# output_dir = '/Users/derek/Dropbox/Master_Bioinformatics/Internship/PANDORA/PANDORA_files/data/outputs/1SJH_1SJE'
    # Change working directory
    os.chdir(output_dir)
    # run Modeller to perform homology modelling
    os.popen('python3 %s > modeller.log' %python_script).read()
    os.chdir(os.path.dirname(PANDORA.PANDORA_path))

    # Parse .log file
    logf = []
    f = open(output_dir + '/modeller.log')
    for line in f:
        if line.startswith(target.id + '.'):
            l = line.split()
            if len(l) > 2:
                logf.append(tuple(l))
    f.close()

    # If keep_IL is true (happens if the target and template are the same), also use the initial model as one of the
    # results. This will also happen while benchmarking.
    if keep_IL:
        # Also take the Initial Loop model. Take the molpdf from the pdb header.
        il_file = [i for i in os.listdir(output_dir) if i.startswith(target.id + '.IL')][0]
        il = open(output_dir + '/' + il_file)
        for line in il:
            if 'MODELLER OBJECTIVE FUNCTION' in line:
                il_molpdf = line.split()[-1]
        f.close()
        # Append the filename and molpdf to the rest of the data
        logf.append((il_file, il_molpdf, ''))

    # Write to output file
    f = open(output_dir + '/molpdf_DOPE.tsv', 'w')
    for i in logf:
        f.write('matched_' + i[0] + ',' + i[1] + ',' + i[2] + '\n')
    f.close()


    # Create Model object of each theoretical model and add it to results
    results = []
    for i in range(len(logf)):
        try:
            m = Model.Model(target, model_path=output_dir + '/' + logf[i][0], output_dir = output_dir,
                                            molpdf=logf[i][1], dope=logf[i][2])

        except:
            print('Something went wrong when calling Model.Model() for case %s' %target.id)
        if benchmark:
            try:
                m.calc_LRMSD(PANDORA.PANDORA_data + '/PDBs/pMHC' + target.MHC_class + '/' + target.id + '.pdb')
                print('l-RMSD for %s: %f' %(target.id, m.lrmsd))
            except:
                print('Something went wrong when calculating l-RMSD for case %s' %target.id)
                pass
            if target.MHC_class == 'II': #only calculate the core L-rmsd for MHCII cases
                try:
                    m.calc_Core_LRMSD(PANDORA.PANDORA_data + '/PDBs/pMHC' + target.MHC_class + '/' + target.id + '.pdb')
                    print('Core l-RMSD for %s: %f' %(target.id, m.core_lrmsd))
                except:
                    print('Something went wrong when calculating core l-RMSD for case %s' %target.id)
                    pass
        results.append(m)


    # Save results as pickle
    if pickle_out:
        pickle.dump(results, open("%s/results_%s.pkl" %(output_dir, os.path.basename(os.path.normpath(output_dir))), "wb"))

    return results
