from Bio.Align import substitution_matrices
PAM30 = substitution_matrices.load('PAM30')
import os
import PANDORA
import pickle
from PANDORA.PMHC import Model
from Bio import Align




def find_template(target, database, seq_based_templ_selection = False):
    ''' Selects the template structure that is best suited as template for homology modelling of the target

    Args:
        target: Target object
        database: Database object
        seq_based_templ_selection: (bool) Use template selection based on template sequences instead of allele.

    Returns: Template object

    '''

    # Sequence based template search if the sequences of the target are provided
    if target.M_chain_seq != '' and seq_based_templ_selection:

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

            return database.MHCI_data[best_template]

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

            return database.MHCII_data[best_template]

    ## For MHC I
    if target.MHC_class == 'I':

        # Find template structures with matching alleles
        putative_templates = {}
        for id in database.MHCI_data:
            if any(x in database.MHCI_data[id].allele_type for x in target.allele_type):
                putative_templates[id] = list(
                    set(target.allele_type) & set(database.MHCI_data[id].allele_type))  # update dict with ID:all matching alleles

        # If the target template already occured in the database, remove it from the dict of putative templates
        putative_templates.pop(target.id)

        # Find the putative template with the best matching peptide
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
        return database.MHCI_data[template_id]


    ## For MHC II
    if target.MHC_class == 'II':

        # Find template structures with matching alleles
        putative_templates = {}
        for id in database.MHCII_data:
            if any(x in database.MHCII_data[id].allele_type for x in target.allele_type):
                # putative_templates[id] = db.MHCII_data[id].allele_type
                putative_templates[id] = list(set(target.allele_type) & set(database.MHCII_data[id].allele_type)) #update dict with ID:all matching alleles

        # If the target template already occured in the database, remove it from the dict of putative templates
        putative_templates.pop(target.id)

        # check if there are any template that have two matching alleles
        # max([len(v) for k,v in putative_templates.items()])

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
        return database.MHCII_data[template_id]

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



def write_modeller_script(target, template, alignment_file, output_dir, n_models=20, stdev=0.1):
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


def run_modeller(output_dir, target, python_script = 'cmd_modeller.py', benchmark = False, pickle_out = True):
    ''' Perform the homology modelling.

    Args:
        output_dir: (string) path to output directory
        target: Target object
        python_script:  (string) path to script that performs the modeller modelling. cmd_modeller.py
        benchmark: (bool) Perform L-RMSD calculations? only works if the target id is an existing pdb id
        pickle_out: (bool) Save a .pkl with the results

    Returns: (list) of Model objects

    '''

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

    # Create Model object of each theoretical model and add it to results
    results = []
    for i in range(len(logf)):
        try:
            m = Model.Model(target, model_path=output_dir + '/' + logf[i][0], output_dir = output_dir,
                                            molpdf=logf[i][1], dope=logf[i][2])
            if benchmark:
                m.calc_LRMSD(PANDORA.PANDORA_data + '/PDBs/pMHC' + target.MHC_class + '/' + target.id + '.pdb')
                m.calc_Core_LRMSD(PANDORA.PANDORA_data + '/PDBs/pMHC' + target.MHC_class + '/' + target.id + '.pdb')
            results.append(m)
        except:
            pass

    # Save results as pickle
    if pickle_out:
        pickle.dump(results, open("%s/results_%s.pkl" %(output_dir, os.path.basename(os.path.normpath(output_dir))), "wb"))

    return results