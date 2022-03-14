from Bio.Align import substitution_matrices
import os
import PANDORA
import pickle
from PANDORA.PMHC import Model
# from Bio import Align
from Bio import pairwise2
from PANDORA.Pandora import Align
import statistics
from Bio.Align import PairwiseAligner
from datetime import datetime


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


def predict_anchors_netMHCIIpan(peptide, allele_type, verbose=True):
    '''Uses netMHCIIpan to predict the binding core of a peptide and infer the anchor positions from that.

    Args:
        target: (Target): Target object containing the peptide sequence and allele type

    Returns: (lst): list of predicted anchor predictions

    '''
    all_netMHCpan_alleles = []
    with open(PANDORA.PANDORA_path + '/../netMHCIIpan-4.0/data/allele.list') as f:
        for line in f:
            all_netMHCpan_alleles.append(line.replace('\n', ''))

    # Format the alles to netMHCIIpan readable format
    target_alleles = [i.split('-')[-1].replace('*', '_') for i in allele_type]

    # The DQ and DP alleles only function in pairs in netMHCIIpan, which we cannot match from our alleles
    # So take the first 3 partially matched allele combinations
    for i in target_alleles:
        if 'DRB' not in i:
            target_alleles = target_alleles + [al for al in all_netMHCpan_alleles if i.replace('_', '') in al][:3]

    # for the DQ and DP cases, alleles are matched (e.g. 'HLA-DQA10102-DQB10602')
    # If two alleles are present is such a combi case, select that combi case as the target allele
    target_alleles_matched = []
    for al in target_alleles:
        hits = 0
        for part in [i.split('*')[-1] for i in allele_type]:
            if part in al:
                hits +=1
        if hits == 2:
            target_alleles_matched.append(al)
    if len(target_alleles_matched) > 0:
        target_alleles = target_alleles_matched

    target_alleles = [i for i in target_alleles if i in all_netMHCpan_alleles]
    # If there are no target alleles that occur in netMHCIIpan, but there is a mouse allele, use all mouse alleles
    # that are supported by netMHCIIpan
    if target_alleles == [] and any(al.startswith('H2') for al in allele_type):
        target_alleles = [i for i in all_netMHCpan_alleles if i.startswith('H-')]

    # If there is no target allele that occurs in netMHCIIpan, just use the standard DRB1_0101
    if target_alleles == []:
        target_alleles = ['DRB1_0101']

    target_alleles_str = ','.join(target_alleles)

    # Setup files
    netmhciipan = PANDORA.PANDORA_path + '/../netMHCIIpan-4.0/netMHCIIpan'
    infile = PANDORA.PANDORA_path + '/../netMHCIIpan-4.0/tmp/%s_%s_%s.txt' %(
        peptide, target_alleles[0], datetime.today().strftime('%Y%m%d_%H%M%S'))
    outfile = PANDORA.PANDORA_path + '/../netMHCIIpan-4.0/tmp/%s_%s_%s_prediction.txt' %(
        peptide, target_alleles[0], datetime.today().strftime('%Y%m%d_%H%M%S'))

    # Write peptide sequence to input file for netMHCIIpan
    with open(infile, 'w') as f:
        f.write(peptide)

    try:
        # run netMHCIIpan
        os.system('%s -f %s -inptype 1 -a %s > %s' % (netmhciipan, infile, target_alleles_str, outfile))

        # Get the output from the netMHCIIpan prediction
        # {allele: (offset, core, core_reliability, score_EL, %rank_EL)}
        pred = {}
        with open(outfile) as f:
            for line in f:
                if peptide in line:
                    ln = [i for i in line[:-1].split(' ') if i != '']
                    pred[ln[1]] = (int(ln[3]), ln[4], float(ln[5]))

        # For each predicted core offset, show the best prediction
        max_scores = [max((i[::-1]) for i in list(pred.values()) if i[0] == s) for s in set([pred[i][0] for i in pred])]
        # order to offset, core, core_reliability
        max_scores = [i[::-1] for i in sorted(max_scores, reverse=True)]

    except ValueError:
        print('Could not predict binding core using netMHCIIpan. Will use the most common anchor positions instead')
        return [3, 6, 8, 11]

    # Remove output file
    os.system('rm %s %s' % (infile, outfile))

    offset, core, core_reliability = max_scores[0]
    # Use the canonical spacing for 9-mer binding cores to predict the anchor positions
    predicted_anchors = [offset + 1, offset + 4, offset + 6, offset + 9]
    # Make sure the prediction is not longer than the peptide just in case
    predicted_anchors = [i for i in predicted_anchors if i <= len(peptide)]

    if verbose:
        print('\tPredicted the binding core using netMHCIIpan (4.0):\n')
        print('\toffset:\t%s\n\tcore:\t%s\n\tprob:\t%s\n' % (offset, core, core_reliability))
        print('\tPredicted peptide anchor residues (assuming canonical spacing): %s' % predicted_anchors)

    return predicted_anchors


def predict_anchors_netMHCpan(peptide, allele_type,
                              verbose=True, rm_output=True):
    '''Uses netMHCIIpan to predict the binding core of a peptide and infer the
    anchor positions from that.

    Args:
        peptide: (str): AA sequence of the peptide
        allele_type: (lst): list of strings of allele types
        verbose: (bool):

    Returns: (lst): list of predicted anchor predictions

    '''
    all_netMHCpan_alleles = []
    with open(PANDORA.PANDORA_path + '/../netMHCpan-4.1/data/allelenames') as f:
        for line in f:
            all_netMHCpan_alleles.append(line.split(' ')[0])#.replace(':',''))
    
    ## Format alleles
    target_alleles = [i.replace('*','') for i in allele_type]
    ## Make sure only netMHCpan available alleles are used
    target_alleles = [i for i in target_alleles if i in all_netMHCpan_alleles]
    
    if len(target_alleles) == 0:
        print('ERROR: The provided Target allele is not available in NetMHCpan-4.1')
        return None
    
    target_alleles_str = ','.join(target_alleles)
    
    # Setup files
    netmhcpan = PANDORA.PANDORA_path + '/../netMHCpan-4.1/netMHCpan'
    infile = PANDORA.PANDORA_path + '/../netMHCpan-4.1/tmp/%s_%s_%s.txt' %(
        peptide, target_alleles[0].replace('*','').replace(':',''), datetime.today().strftime('%Y%m%d_%H%M%S'))
    outfile = PANDORA.PANDORA_path + '/../netMHCpan-4.1/tmp/%s_%s_%s_prediction.txt' %(
        peptide, target_alleles[0].replace(':',''), datetime.today().strftime('%Y%m%d_%H%M%S'))
    
    # Write peptide sequence to input file for netMHCIIpan
    with open(infile, 'w') as f:
        f.write(peptide)
    
    os.system('%s -p %s -a %s > %s' %(netmhcpan, infile, target_alleles_str, outfile))
    
    # Get the output from the netMHCIIpan prediction
    # {allele: (core, %rank_EL)}
    pred = {}
    with open(outfile) as f:
        for line in f:
            if peptide in line and not line.startswith('#'):
                ln = [i for i in line[:-1].split(' ') if i != '']
                #ln[3] is core, ln[9] is Icore
                try:
                    pred[ln[1]].append((ln[3], float(ln[12])))
                except KeyError:
                    pred[ln[1]] = [(ln[3], float(ln[12]))]

    # Sort each allele result per Rank_EL
    for allele in pred:
        pred[allele] = list(sorted(pred[allele], key=lambda x:x[1]))
    
    if len(pred) == 0:
        print('ERROR: NetMHCpan-4.1 was not able to find any binding core for')
        print('the provided peptide and MHC allele')
        return None
    
    # For every allele, the binding core is predicted. Take the allele with the highest reliability score
    best_allele = min((pred[i][0][1], i) for i in pred)[1]
    
    # Do a quick alignment of the predicted core and the peptide to find the anchors. (the predicted binding core can
    # contain dashes -. Aligning them makes sure you take the right residue as anchor.
    alignment = pairwise2.align.globalxx(peptide, pred[best_allele][0][0])
    
    #If there are multiple possible solutions, take the one with no gap at the anchor (second) position
    if len(alignment)>1:
        flag = False
        #Search for the options without gap in the second postions
        for prediction in alignment:
            if prediction[1][1] != '-' and prediction[0][1] != '-':
                pept1 = prediction[0]
                pept2 = prediction[1]
                flag = True
                break
        #If no options are available, take the first one
        if flag==False:
             pept1 = alignment[0][0]
             pept2 = alignment[0][1]
        
    else:
        pept1 = alignment[0][0]
        pept2 = alignment[0][1]
    
    # Remove gaps if in the same position
    to_remove = []
    for i, (aa1, aa2) in enumerate(zip(pept1, pept2)):
        if aa1 == aa2 == '-' and i != 0:
            to_remove.append(i)
    for x in reversed(to_remove):
       pept1 = pept1[0:x:]+pept1[x+1::]
       pept2 = pept2[0:x:]+pept2[x+1::]

    if verbose:
        print('Query peptide aligned to the core:')
        print(pept1)
        print(pept2)
    
    # Find the anchors by finding the first non dash from the left and from the right
    # Define chanonical ancors as starting list
    predicted_anchors = [2,len(peptide)]

    # Find the first anchor
    p1 = 0
    p2 = 0
    for i in range(len(pept2)):
        # if the second position has no gaps
        if i == 1 and pept2[i] != '-' and pept1[i] != '-':
            predicted_anchors[0] = p1 + 1
            break
        elif i > 1 and pept2[i] != '-':
            predicted_anchors[0] = p1 + 1
            break
        if pept1[i] != '-':
            p1 += 1
        if pept2[i] != '-':
            p2 += 1
    
    # Find the second anchor
    for i in range(len(pept2)):
        if pept2[::-1][i] != '-':
            predicted_anchors[1] = len([j for j in pept1[:len(pept1) -i] if j != '-'])
            #predicted_anchors[1] = len([j for j in pept2[::-1][i] if j != '-'])
            break
    
    if verbose:
        print('\tPredicted the binding core using netMHCpan (4.1):\n')
        print('\tIcore:\t%s\n\t%%Rank EL:\t%s\n' %(pred[best_allele][0][0], pred[best_allele][0][1] ))
        print('\tPredicted peptide anchor residues (assuming canonical spacing): %s' %predicted_anchors)
    
    if rm_output:
        os.system('rm %s' %infile)
        os.system('rm %s' %outfile)

    return predicted_anchors


def score_peptide_alignment_MHCI(target, template, substitution_matrix='PAM30'):
    ''' Calculate the alignment score of the target and template peptide

    Args:
        target: (Target): Target object
        template: (Template): Template object
        substitution_matrix: (str): name of subtitution matrix, default is PAM30 (BLOSUM80 etc)

    Returns: (flt): alignment score

    '''
    # Dario don't worry, I didn't change the code, I just moved it to a function, so peptide similarity can be
    # calculated for user defined templates as well.
    substitution_matrix = substitution_matrices.load(substitution_matrix)
    score = 0
    try:
        pept_anchs = target.anchors
    except:
        pept_anchs = [1, len(target.peptide) - 1]

    temp_pept = template.peptide
    temp_anchs = template.anchors
    aligned_pept, aligned_temp_pept = align_peptides(target.peptide,
                                                     pept_anchs[0], pept_anchs[1],
                                                     temp_pept,
                                                     temp_anchs[0], temp_anchs[1])

    aligned_pept = aligned_pept.replace('-', '*')
    aligned_temp_pept = aligned_temp_pept.replace('-', '*')
    # min_len = min([len(target.peptide), len(temp_pept)])
    # score -= ((abs(len(target.peptide) - len(temp_pept)) ** 2.4))  # !!!  ## Gap Penalty #Is now handled by normal PAM30
    for i, (aa, bb) in enumerate(zip(aligned_pept, aligned_temp_pept)):
        try:
            # gain = MatrixInfo.pam30[aa, bb]
            gain = substitution_matrix[aa, bb]
            score += gain
        except KeyError:
            try:
                # gain = MatrixInfo.pam30[bb, aa]
                gain = substitution_matrix[bb, aa]
                score += gain
            except KeyError:
                score = -50
                pass

    return score


def score_peptide_alignment_MHCII(target, template, substitution_matrix='PAM30'):
    ''' Calculate the alignment score of the target and template peptide using pairwise alignment

    Args:
        target: (Target): Target object
        template: (Template): Template object
        substitution_matrix: (str): name of subtitution matrix, default is PAM30 (BLOSUM80 etc)

    Returns: (flt): alignment score

    '''

    # define the peptide and p1 anchor position
    temp_pept = template.peptide
    temp_p1 = template.anchors[0]
    tar_pept = target.peptide
    tar_p1 = target.anchors[0]

    # align based on first anchor position, fill in the ends with '-' to make them equal length
    temp_pept = '*' * (tar_p1 - temp_p1) + temp_pept
    tar_pept = '*' * (temp_p1 - tar_p1) + tar_pept
    temp_pept = temp_pept + '*' * (len(tar_pept) - len(temp_pept))
    tar_pept = tar_pept + '*' * (len(temp_pept) - len(tar_pept))

    # Perform a pairwise alignment. Make sure no gaps are introduced.
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load(substitution_matrix)
    aligner.gap_score = -1000
    aligner.end_open_gap_score = -1000000
    aligner.internal_open_gap_score = -10000

    # Align the sequences
    aligned = aligner.align(tar_pept, temp_pept)

    return aligned.score


def find_template(target, database, best_n_templates = 1, benchmark=False):
    ''' Selects the template structure that is best suited as template for homology modelling of the target

    Args:
        target: Target object
        database: Database object
        seq_based_templ_selection: (bool) Use template selection based on template sequences instead of allele.

    Returns: Template object

    '''

    ## For MHC I
    if target.MHC_class == 'I':

        # Define available alleles in database
        available_alleles = []
        for ID in database.MHCI_data:
            if benchmark and ID == target.id:
                pass
            else:
                available_alleles.extend(database.MHCI_data[ID].allele_type)
        available_alleles = list(set(available_alleles))

        # Adapt the target allele name if necessary
        #target_alleles = [allele_name_adapter(allele, available_alleles) for allele in target.allele_type]
        target_alleles = allele_name_adapter(target.allele_type, available_alleles)
        target_alleles = list(set(target_alleles))

        # Find template structures with matching alleles
        putative_templates = {}
        for ID in database.MHCI_data:
            if benchmark and ID == target.id:
                pass
            else:
                for tar_allele in target_alleles:
                    if any(tar_allele in put_temp_allele for put_temp_allele in database.MHCI_data[ID].allele_type):
                        # update dict with ID:all matching alleles
                        #TODO: is this list of matching allele obsolete?
                        putative_templates[ID] = list(
                            set(target.allele_type) & set(database.MHCI_data[ID].allele_type))

        # If the target template already occured in the database, remove it from the dict of putative templates
        #putative_templates.pop(target.id)

        # Find the putative template with the best matching peptide
        pos_list = []
        for ID in putative_templates:
            score = score_peptide_alignment_MHCI(target, database.MHCI_data[ID], substitution_matrix='PAM30')
            pos_list.append((score, database.MHCI_data[ID].peptide, ID))

        if len(pos_list) == 0:
            raise Exception('Pandora could not find any putative template! Please try to define your own template or contact us for help')
        # Take the putative template with the max scoring peptide
        # template_id = pos_list[[i[0] for i in pos_list].index(max([i[0] for i in pos_list]))][2]
        # Return the Template object of the selected template that will be used for homology modelling

        template_id = [i[-1] for i in sorted(pos_list, key=lambda elem: elem[0], reverse=True)][:best_n_templates]
        scores = sorted(pos_list, key=lambda elem: elem[0], reverse=True)[:best_n_templates]

        templates = [database.MHCI_data[tmpl] for tmpl in template_id]
        keep_IL = any(check_target_template(target, tmpl) for tmpl in templates)

        return templates, scores, keep_IL


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

        # Make sure there is no template with only 3 anchors for benchmarking.
        if benchmark:
            putative_templates = {k:v for k,v in putative_templates.items() if len(database.MHCII_data[k].anchors) == 4}


        # Find the peptide with the highest alignment score. If there are multiple, take the first one with same
        # same anchor positions
        # template_id = find_best_template_peptide(target=target,
        #                                          templates=[database.MHCII_data[i] for i in putative_templates])

        # Find the putative template with the best matching peptide
        pos_list = []
        for ID in putative_templates:
            score = score_peptide_alignment_MHCII(target, database.MHCII_data[ID], substitution_matrix='PAM30')
            pos_list.append((score, database.MHCII_data[ID].peptide, ID))

        if len(pos_list) == 0:
            raise Exception('Pandora could not find any putative template! Please try to define your own template or contact us for help')
        # Take the putative template with the max scoring peptide
        # template_id = pos_list[[i[0] for i in pos_list].index(max([i[0] for i in pos_list]))][2]        # Return the Template object of the selected template that will be used for homology modelling


        template_id = [i[-1] for i in sorted(pos_list, key=lambda elem: elem[0], reverse=True)][:best_n_templates]
        scores = sorted(pos_list, key=lambda elem: elem[0], reverse=True)[:best_n_templates]

        templates = [database.MHCII_data[tmpl] for tmpl in template_id]
        keep_IL = any(check_target_template(target, tmpl) for tmpl in templates)

        return templates, scores, keep_IL

            # return database.MHCII_data[template_id], check_target_template(target, database.MHCII_data[template_id])

    # # Sequence based template search if the sequences of the target are provided
    # elif target.M_chain_seq != '' and seq_based_templ_selection:
    #
    #     if target.MHC_class == 'I':
    #
    #         # define target sequences
    #         tar_seq = database.MHCI_data[target.id].M_chain_seq
    #         tar_pept = database.MHCI_data[target.id].peptide
    #         # keep track of alignment scores
    #         scores = {}
    #         # Perform a pairwise alignment of the target and all templates for the MHC M chain and peptide
    #         for i in database.MHCI_data:
    #             aligner = Align.PairwiseAligner()
    #             aligner.substitution_matrix = substitution_matrices.load("BLOSUM80")  # PAM30 for pept??
    #
    #             M_score = aligner.align(tar_seq, database.MHCI_data[i].M_chain_seq).score
    #             P_score = aligner.align(tar_pept, database.MHCI_data[i].peptide).score
    #
    #             scores[i] = (M_score, P_score)
    #         # Remove the target structure from this dict, you cannot select the target as template
    #         scores.pop(target.id, None)
    #         # take the 10 best scoring templates
    #         best_MHCs = sorted(scores, key=scores.get, reverse=True)[:10]
    #         # take the template with the best scoring peptide
    #         best_template = max((v[1], k) for k, v in scores.items() if k in best_MHCs)[1]
    #
    #         return database.MHCI_data[best_template], check_target_template(target, database.MHCI_data[best_template])
    #
    #     if target.MHC_class == 'II':
    #         # define target sequences
    #         tar_seq = database.MHCII_data[target.id].M_chain_seq + database.MHCII_data[target.id].N_chain_seq
    #         tar_pept = database.MHCII_data[target.id].peptide
    #         # keep track of alignment scores
    #         scores = {}
    #
    #         for i in database.MHCII_data:
    #             aligner = Align.PairwiseAligner()
    #             aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")  # or PAM30 ??
    #
    #             temp_seq = database.MHCII_data[i].M_chain_seq + database.MHCII_data[i].N_chain_seq
    #             MN_score = aligner.align(tar_seq, temp_seq).score
    #             P_score = aligner.align(tar_pept, database.MHCII_data[i].peptide).score
    #
    #             scores[i] = (MN_score, P_score)
    #         # Remove the target structure from this dict, you cannot select the target as template
    #         scores.pop(target.id, None)
    #         # take the 10 best scoring templates
    #         best_MHCs = sorted(scores, key=scores.get, reverse=True)[:10]
    #         # take the template with the best scoring peptide
    #         best_template = max((v[1], k) for k, v in scores.items() if k in best_MHCs)[1]
    #
    #         return database.MHCII_data[best_template], check_target_template(target, database.MHCI_data[best_template])


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
                if 'self.residue_range' in line and 'M.selection' in line:
                    myloopscript.write(line % (anch[0]+1, anch[-1]-1))
                elif 'SPECIAL_RESTRAINTS_BREAK' in line:
                    break
                elif 'contact_file = open' in line:
                    myloopscript.write(line %target.id)
                else:
                    myloopscript.write(line)
            MyL_temp.close()

    if target.MHC_class == 'II':
        with open(output_dir + '/MyLoop.py', 'w') as myloopscript:
            MyL_temp = open(PANDORA.PANDORA_path + '/Pandora/MyLoop_template_II.py', 'r')
            for line in MyL_temp:
                if 'self.residue_range' in line and 'M.selection' in line:
                    if anch[0] == 0:
                        anch_1 = 1
                    else:
                        anch_1 = anch[0]
                    if anch[-1] == (len(target.peptide)-1):
                        anch_term = len(target.peptide)
                    else:
                        anch_term = anch[-1]
                    #Write first and last anchors, to keep only the flanking regions flexible
                    myloopscript.write(line % (1, anch_1, anch_term, len(target.peptide)))
                    #for i in range(len(anch)-1): # Write all the inbetween acnhors if they are there
                    #    myloopscript.write(line % (anch[i] + 2, anch[i+1]))
                    #myloopscript.write(line % (anch[-1] + 2, len(target.peptide))) # Write the last anchor
                elif 'SPECIAL_RESTRAINTS_BREAK' in line:
                    break
                elif 'contact_file = open' in line:
                    myloopscript.write(line %target.id)
                else:
                    myloopscript.write(line)
            MyL_temp.close()

    with open(output_dir.replace('\\ ', ' ') + '/cmd_modeller_ini.py', 'w') as modscript:
        cmd_m_temp = open(PANDORA.PANDORA_path + '/Pandora/cmd_modeller_ini.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line % os.path.basename(alignment_file))
            elif 'knowns' in line:
                if type(template)==list:
                    modscript.write(
                        'knowns = (%s), sequence = "%s",\n' % (','.join(['"' + i.id + '"' for i in template]), target.id))
                else:
                    modscript.write(
                        'knowns = (%s), sequence = "%s",\n' % ('"' + template.id + '"', target.id))
                # modscript.write(line % ('(' + ','.join([i.id for i in template]) + ')', target.id))
            else:
                modscript.write(line)
        cmd_m_temp.close()

# alignment_file = mod.alignment.alignment_file
# output_dir = mod.output_dir
# template = mod.template
# helix = [3, 8]
# BETA-SHEET-MARKER


def write_modeller_script(target, template, alignment_file, output_dir, n_homology_models=1, n_loop_models = 20,
                          loop_refinement='slow', n_jobs=None, stdev=0.1, helix = False, sheet = False):
    ''' Write script that refines the loops of the peptide
    
    Args:
        target (PANDORA.PMHC.PMHC.Target): Target object
        template (PANDORA.PMHC.PMHC.Template): Template object
        alignment_file (str): path to alignment file
        output_dir (str): path to output directory
        n_homology_models (int): number of homology models that are generated per run.
        n_loop_models (int): number of loop models modeller generates per homology model
        n_jobs (int): number of parallel jobs. Is recommended to use at most as many jobs as the number of models:
            ore will not add any benefit but might occupy cores unnecessarily.
        loop_refinement (str): Level of loop refinement: very_fast,fast,slow,very_slow,slow_large.
            Defaults to slow
        stdev (float): standard deviation of modelling restraints. Higher = more flexible restraints.
        helix (list): List of the alpha helix start and end-positions as integers. I.e. [3,8] for a helix between
            peptide residue 3 and 8.
        sheet (list): List containing: start position of B-sheet 1, start position of B-sheet 2 and the length of the
            B-sheet in h-bonds. For example: ["O:2:P","N:54:M",2] for a parallel B-sheet; The sheet starts
            at the Oxigen atom of the 2nd residue of chain P and at the Nitrogen of the 54th residue of
            chain M and has a length of 2 H-bonds. Or; ["N:6:P", "O:13:P", -3], with -3 denoting an
            anti-parallel B-sheet with a length of 3 H-bonds.

    '''

    anch = target.anchors

    if target.MHC_class == 'I':
        with open(output_dir.replace('\\ ', ' ') + '/MyLoop.py', 'w') as myloopscript:
            MyL_temp = open(PANDORA.PANDORA_path + '/Pandora/MyLoop_template.py', 'r')
            for line in MyL_temp:
                if 'self.residue_range' in line and 'M.selection' in line:
                    myloopscript.write(line %(anch[0]+1, anch[-1]-1))  # write the first anchor
                elif 'contact_file = open' in line:
                    myloopscript.write(line %(target.id))
                elif 'STDEV MARKER' in line:
                    myloopscript.write(line %(stdev))
                elif helix and 'ALPHA-HELIX-MARKER' in line:
                    myloopscript.write(line.replace('# ALPHA-HELIX-MARKER', 'rsr.add(M.secondary_structure.alpha(self.residue_range("%s:P", "%s:P")))' %(helix[0], helix[1])))
                elif sheet and 'BETA-SHEET-MARKER' in line:
                    myloopscript.write(line.replace('# BETA-SHEET-MARKER', 'rsr.add(M.secondary_structure.sheet(atoms["%s"], atoms["%s"], sheet_h_bonds=%s))' %(sheet[0], sheet[1], sheet[2])))
                else:
                    myloopscript.write(line)
            MyL_temp.close()

    if target.MHC_class == 'II':
        with open(output_dir.replace('\\ ', ' ') + '/MyLoop.py', 'w') as myloopscript:
            MyL_temp = open(PANDORA.PANDORA_path + '/Pandora/MyLoop_template_II.py', 'r')
            for line in MyL_temp:
                if 'self.residue_range' in line and 'M.selection' in line:
                    if anch[0] == 0:
                        anch_1 = 1
                    else:
                        anch_1 = anch[0]
                    if anch[-1] == (len(target.peptide)-1):
                        anch_term = len(target.peptide)
                    else:
                        anch_term = anch[-1]
                    #Write first and last anchors, to keep only the flanking regions flexible
                    myloopscript.write(line % (1, anch_1, anch_term, len(target.peptide)))
                    #for i in range(len(anch)-1): # Write all the inbetween acnhors if they are there
                    #    myloopscript.write(line % (anch[i] + 2, anch[i+1]))
                    #myloopscript.write(line % (anch[-1] + 2, len(target.peptide))) # Write the last anchor
                elif 'contact_file = open' in line:
                    myloopscript.write(line %(target.id))
                elif 'STDEV MARKER' in line:
                    myloopscript.write(line %(stdev))
                elif helix and 'ALPHA-HELIX-MARKER' in line:
                    myloopscript.write(line.replace('# ALPHA-HELIX-MARKER', 'rsr.add(M.secondary_structure.alpha(self.residue_range("%s:P", "%s:P")))' %(helix[0], helix[1])))
                elif sheet and 'BETA-SHEET-MARKER' in line:
                    myloopscript.write(line.replace('# BETA-SHEET-MARKER', 'rsr.add(M.secondary_structure.sheet(atoms["%s"], atoms["%s"], sheet_h_bonds=%s))' %(sheet[0], sheet[1], sheet[2])))
                else:
                    myloopscript.write(line)
            MyL_temp.close()

    with open(output_dir.replace('\\ ', ' ') + '/cmd_modeller.py', 'w') as modscript:
        cmd_m_temp = open(PANDORA.PANDORA_path + '/Pandora/cmd_modeller_template.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line %(os.path.basename(alignment_file)))
            elif 'knowns' in line:
                if type(template)==list:
                    modscript.write(
                        'knowns = (%s), sequence = "%s",\n' % (','.join(['"' + i.id + '"' for i in template]), target.id))
                else:
                    modscript.write(
                        'knowns = (%s), sequence = "%s",\n' % ('"' + template.id + '"', target.id))
                # modscript.write(line %(','.join([i.id for i in template]), target.id))
            elif 'a.ending_model' in line:
                modscript.write(line % (n_homology_models))
            elif 'a.loop.ending_model' in line:
                modscript.write(line % (n_loop_models))
            elif 'a.loop.md_level' in line:
                modscript.write('a.loop.md_level = MA.refine.%s   # Loop model refinement level' %(loop_refinement))
            else:
                if n_jobs != None: #If this is a parallel job
                    if 'PARALLEL_JOB_LINE_TO_COMPLETE' in line:
                        modscript.write(line %(str(n_jobs))) #specify the number of cores
                    else:
                        modscript.write(line)  #Write the line as it is
                else: #If this is not a parallel job
                    if 'PARALLEL_JOB_LINE' in line: #do not write the lines requested for parallelization
                        pass
                    else:
                        modscript.write(line)  #Write the line as it is
        cmd_m_temp.close()



def run_modeller(output_dir, target, python_script = 'cmd_modeller.py', benchmark = False, pickle_out = True,
                 keep_IL = False, RMSD_atoms = ['C', 'CA', 'N', 'O']):
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
        if line.startswith(target.id + '.'):   #target.id
            l = line.split()
            if len(l) == 3: #TODO: make sure the line is reporting the model with tis score. Format: model, molpdf, dope.
                logf.append(tuple(l))
    f.close()

    # If keep_IL is true (happens if the target and template are the same), also use the initial model as one of the
    # results. This will also happen while benchmarking.
    if keep_IL:
        # Also take the Initial Loop model. Take the molpdf from the pdb header.
        il_file = [i for i in os.listdir(output_dir) if i.startswith(target.id + '.IL')][0]
        # il = open(output_dir + '/' + il_file)
        # for line in il:
        #     if 'MODELLER OBJECTIVE FUNCTION' in line:
        #         il_molpdf = line.split()[-1]
        # f.close()
        # Create a fake molpdf/dope score for the IL model: the best molpdf/dope from the real models - 1
        try:
            fake_molpdf = str(min(float(i[1]) for i in logf) - 1)
            fake_dope = str(min(float(i[2]) for i in logf) - 1)
        except ValueError:
            fake_molpdf = -10000
            fake_dope = -10000
            print('WARNING: ValueError exception raised while assigning fake molpdf and dope to IL model')
        # Append the filename and molpdf to the rest of the data
        logf.append((il_file, fake_molpdf, fake_dope))

    # Sort output by molpdf
    logf.sort(key=lambda tup:float(tup[1]))

    # Write to output file
    f = open(output_dir + '/molpdf_DOPE.tsv', 'w')
    for i in logf:
        f.write(i[0] + '\t' + i[1] + '\t' + i[2] + '\n')
    f.close()


    # Create Model object of each theoretical model and add it to results
    results = []
    for i in range(len(logf)):
        try:
            m = Model.Model(target, model_path=output_dir + '/' + logf[i][0], output_dir = output_dir,
                                            molpdf=logf[i][1], dope=logf[i][2])

        except:
            print('WARNING: Error raised while calling Model.Model() for case %s' %target.id)
            
        # if benchmark:
        #     try:
        #         m.calc_LRMSD(PANDORA.PANDORA_data + '/PDBs/pMHC' + target.MHC_class + '/' + target.id + '.pdb',
        #                      atoms = RMSD_atoms)
        #         # print('l-RMSD for %s: %f' %(target.id, m.lrmsd))
        #     except:
        #         print('Something went wrong when calculating l-RMSD for case %s' %target.id)
        #         pass
        #     if target.MHC_class == 'II': #only calculate the core L-rmsd for MHCII cases
        #         try:
        #             m.calc_Core_LRMSD(PANDORA.PANDORA_data + '/PDBs/pMHC' + target.MHC_class + '/' + target.id + '.pdb',
        #                      atoms = RMSD_atoms)
        #             # print('Core l-RMSD for %s: %f' %(target.id, m.core_lrmsd))
        #         except:
        #             print('Something went wrong when calculating core l-RMSD for case %s' %target.id)
        #             pass
        results.append(m)


    # Save results as pickle
    if pickle_out:
        pickle.dump(results, open("%s/results_%s.pkl" %(output_dir, os.path.basename(os.path.normpath(output_dir))), "wb"))

    return results

def align_peptides(seq1, anch1_seq1, anch2_seq1, seq2, anch1_seq2, anch2_seq2):
    '''
    Align two MHC-I peptides making overlap the anchors.
    This function does NOT use an alignment matrix (e.g. BLOSUM, PAM, etc).
    It computes a simple anchor position alignment and inserts gap in the
    middle part to make the final sequences have the same lenghts.

    Args:
        seq1(str) : sequence of the first peptide.
        anch1_seq1(int) : position of the first anchor of seq1. Position must be given in Python numbering (0-N)
        anch2_seq1(int) : position of the second anchor of seq1. Position must be given in Python numbering (0-N)
        seq2(str) : sequence of the second peptide.
        anch1_seq1(int) : position of the first anchor of seq1. Position must be given in Python numbering (0-N)
        anch2_seq1(int) : position of the second anchor of seq1. Position must be given in Python numbering (0-N)

    Returns:
        ali_seq1(str)
    '''

    seq1_core = anch2_seq1 - anch1_seq1
    seq2_core = anch2_seq2 - anch1_seq2
    tail1 = [x for x in seq1[anch2_seq1:]]
    tail2 = [x for x in seq2[anch1_seq2:]]

    list1 = [x for x in seq1]
    list2 = [x for x in seq2]
    #Adding gaps in cores
    if seq1_core > seq2_core:
        for x in range(seq1_core - seq2_core):
            list2.insert(int(len(seq2)/2), '-')
    elif seq1_core < seq2_core:
        for x in range(seq2_core - seq1_core):
            list1.insert(int(len(seq1)/2), '-')
    ### Adding gaps in heads
    if anch1_seq1 > anch1_seq2:
        for x in range(anch1_seq1 - anch1_seq2):
            list2.insert(0, '-')
    elif anch1_seq1 < anch1_seq2:
        for x in range(anch1_seq2 - anch1_seq1):
            list1.insert(0, '-')
    ### Adding gaps in heads
    if len(tail1) > len(tail2):
        for x in range(len(tail1) - len(tail2)):
            list2.insert(-1, '-')
    elif len(tail1) < len(tail2):
        for x in range(len(tail1) - len(tail2)):
            list1.insert(-1, '-')

    ali_seq1 = ('').join(list1)
    ali_seq2 = ('').join(list2)
    return ali_seq1, ali_seq2

def allele_name_adapter(allele, available_alleles):
    '''
    Cuts the given allele name to make it consistent with the alleles in allele_ID.

    Args:
        allele(list) : Allele names
        allele_ID(dict) : Dictionary of structure IDs (values) in the dataset for each allele (keys)
        
    Returns:
        allele(list) : List of adapted (cut) allele names
    '''
    #homolog_allele = '--NONE--'
    for a in range(len(allele)):
        if allele[a].startswith('HLA'):      # Human
            if any(allele[a] in key for key in list(available_alleles)):
                pass
            elif any(allele[a][:8] in key for key in list(available_alleles)):
                allele[a] = allele[a][:8]
            elif any(allele[a][:6] in key for key in list(available_alleles)):
                allele[a] = allele[a][:6]
            else:
                allele[a] = allele[a][:4]
        elif allele[a].startswith('H2'):    # Mouse
            #homolog_allele = 'RT1'
            if any(allele[a] in key for key in list(available_alleles)):
                pass
            elif any(allele[a][:4] in key for key in list(available_alleles)):
                allele[a] = allele[a][:4]
            else:
                allele[a] = allele[a][:3]
        elif allele[a].startswith('RT1'):          # Rat
            #homolog_allele = 'H2'
            if any(allele[a] in key for key in list(available_alleles)):
                pass
            elif any(allele[a][:5] in key for key in list(available_alleles)):
                allele[a] = allele[a][:5]
            else:
                allele[a] = allele[a][:4]
        elif allele[a].startswith('BoLA'):        # Bovine
            if any(allele[a] in key for key in list(available_alleles)):
                pass
            elif any(allele[a][:10] in key for key in list(available_alleles)):
                allele[a] = allele[a][:10]
            elif any(allele[a][:7] in key for key in list(available_alleles)):
                allele[a] = allele[a][:7]
            else:
                allele[a] = allele[a][:5]
        elif allele[a].startswith('SLA'):        # Suine
            if any(allele[a] in key for key in list(available_alleles)):
                pass
            elif any(allele[a][:9] in key for key in list(available_alleles)):
                allele[a] = allele[a][:9]
            elif any(allele[a][:6] in key for key in list(available_alleles)):
                allele[a] = allele[a][:6]
            else:
                allele[a] = allele[a][:4]
        elif allele[a].startswith('MH1-B'):        # Chicken
            if any(allele[a] in key for key in list(available_alleles)):
                pass
            elif any(allele[a][:8] in key for key in list(available_alleles)):
                allele[a] = allele[a][:8]
            else:
                allele[a] = allele[a][:6]
        elif allele[a].startswith('MH1-N'):        # Chicken
            if any(allele[a] in key for key in list(available_alleles)):
                pass
            elif any(allele[a][:9] in key for key in list(available_alleles)):
                allele[a] = allele[a][:9]
            else:
                allele[a] = allele[a][:6]
        elif allele[a].startswith('BF2'):        # Chicken
            if any(allele[a] in key for key in list(available_alleles)):
                pass
            elif any(allele[a][:6] in key for key in list(available_alleles)):
                allele[a] = allele[a][:6]
            else:
                allele[a] = allele[a][:4]
        elif allele[a].startswith('Mamu'):       # Monkey
            if any(allele[a] in key for key in list(available_alleles)):
                pass
            elif any(allele[a][:13] in key for key in list(available_alleles)):
                allele[a] = allele[a][:13]
            elif any(allele[a][:9] in key for key in list(available_alleles)):
                allele[a] = allele[a][:9]
            else:
                allele[a] = allele[a][:5]
        elif allele[a].startswith('Eqca'):        # Horse
            if any(allele[a] in key for key in list(available_alleles)):
                pass
            elif any(allele[a][:10] in key for key in list(available_alleles)):
                allele[a] = allele[a][:10]
            elif any(allele[a][:7] in key for key in list(available_alleles)):
                allele[a] = allele[a][:7]
            else:
                allele[a] = allele[a][:5]
    return(allele)#, homolog_allele)
