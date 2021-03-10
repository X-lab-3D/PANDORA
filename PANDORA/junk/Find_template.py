
from Bio.Align import substitution_matrices
PAM30 = substitution_matrices.load('PAM30')


def find_template(target, database):
    ''' Selects the structure that is best suited as template for homology modelling of the target

    :param target: (Target) object
    :param database: (Database) object
    :return: (Template) Template object of the best structure
    '''


    ## For MHC I
    if target.MHC_class == 'I':

        # Find template structures with matching alleles
        putative_templates = {}
        for id in database.MHCI_data:
            if any(x in database.MHCI_data[id].allele for x in target.allele):
                putative_templates[id] = list(
                    set(target.allele) & set(database.MHCI_data[id].allele))  # update dict with ID:all matching alleles

        # If the target template already occured in the database, remove it from the dict of putative templates
        putative_templates.pop(target.PDB_id)

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
            if any(x in database.MHCII_data[id].allele for x in target.allele):
                # putative_templates[id] = db.MHCII_data[id].allele
                putative_templates[id] = list(set(target.allele) & set(database.MHCII_data[id].allele)) #update dict with ID:all matching alleles

        # If the target template already occured in the database, remove it from the dict of putative templates
        putative_templates.pop(target.PDB_id)

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






#
# def select_template(IDD, target_id, allele, length, pept, print_results): #homolog_allele,
#     '''
#     Select template for p:MHC I Modelling.
#
#     Args:
#         IDD(dict) : dictionary containing all the templates ID and relative informations
#         allele(str) : target allele name
#         length(int) : target peptide length
#         pept(str) : target peptide sequence
#         print_result(bool) : if True prints out the selected template structure and peptide
#
#     '''
#
#     putative_templates = []
#     max_pos = -1000
#     pos_list = []
#
#     for ID in IDD:
#         for a in allele:
#             if any(a in key for key in IDD[ID]['allele']): #or homolog_allele in IDD[ID]['allele']:                       ## Same Allele
#                 putative_templates.append(ID)
#     putative_templates = list(set(putative_templates))
#
#
#     for ID in putative_templates:
#         score = 0
#         temp_pept = IDD[ID]['pept_seq']
#         min_len = min([length, len(temp_pept)])
#         score -= ((abs(length - len(temp_pept)) ** 2.4)) #!!!  ## Gap Penalty
#         for i, (aa, bb) in enumerate(zip(pept[:min_len], temp_pept[:min_len])):
#             try:
#                 #gain = MatrixInfo.pam30[aa, bb]
#                 gain = PAM30[aa, bb]
#                 score += gain
#             except KeyError:
#                 try:
#                     #gain = MatrixInfo.pam30[bb, aa]
#                     gain = PAM30[bb, aa]
#                     score += gain
#                 except KeyError:
#                     score = -50
#                     pass
#
#         if score > max_pos:
#             max_pos = score
#         pos_list.append((score, temp_pept, ID))
#
#     max_list = []
#     for pos in pos_list:
#         if pos[0] == max_pos:
#             max_list.append(pos)
#     #maxs.append(max_pos)
#     #maxsl.append(len(max_list))
#
#     if len(max_list) == 0:
#         return (target_id, pept, allele, "NA", "NA", "NA", "No positive scoring template peptides")
#     elif len(max_list) == 1:
#         template = max_list[0]
#         template_ID = template[2]
#         template_pept = template[1]
#     else:
#         template = (choice(max_list))
#         template_ID = template[2]
#         template_pept = template[1]
#     if print_results:
#         print('####################################################')
#         print('')
#         print('Peptide:  ', template[1])
#         print('')
#         print('Pos:  ', template)
#         print('')
#         print('####################################################')
#
#     return((template, template_ID, template_pept))
#
#
# def find_template(Target, Database): #homolog_allele,  #TODO Dario's function
#     ''' Select template for p:MHC I Modelling.
#
#     Args:
#         IDD(dict) : dictionary containing all the templates ID and relative informations
#         allele(str) : target allele name
#         length(int) : target peptide length
#         pept(str) : target peptide sequence
#         print_result(bool) : if True prints out the selected template structure and peptide
#
#     '''
#
#     putative_templates = []
#     max_pos = -1000
#     pos_list = []
#
#     for ID in IDD:
#         for a in allele:
#             if any(a in key for key in IDD[ID]['B_allele']): #or homolog_allele in IDD[ID]['allele']:                       ## Same Allele
#                 putative_templates.append(ID)
#     putative_templates = list(set(putative_templates))
#
#
#     for ID in putative_templates:
#         score = 0
#         temp_pept = IDD[ID]['C']
#         min_len = min([length, len(temp_pept)])
#         score -= ((abs(length - len(temp_pept)) ** 2.4)) #!!!  ## Gap Penalty
#         for i, (aa, bb) in enumerate(zip(pept[:min_len], temp_pept[:min_len])):
#             try:
#                 #gain = MatrixInfo.pam30[aa, bb]
#                 gain = PAM30[aa, bb]
#                 score += gain
#             except KeyError:
#                 try:
#                     #gain = MatrixInfo.pam30[bb, aa]
#                     gain = PAM30[bb, aa]
#                     score += gain
#                 except KeyError:
#                     score = -50
#                     pass
#
#         if score > max_pos:
#             max_pos = score
#         pos_list.append((score, temp_pept, ID))
#
#     max_list = []
#     for pos in pos_list:
#         if pos[0] == max_pos:
#             max_list.append(pos)
#     #maxs.append(max_pos)
#     #maxsl.append(len(max_list))
#
#     if len(max_list) == 0:
#         return (target_id, pept, allele, "NA", "NA", "NA", "No positive scoring template peptides")
#     elif len(max_list) == 1:
#         template = max_list[0]
#         template_ID = template[2]
#         template_pept = template[1]
#     else:
#         template = (choice(max_list))
#         template_ID = template[2]
#         template_pept = template[1]
#     # if print_results:
#     #     print('####################################################')
#     #     print('')
#     #     print('Peptide:  ', template[1])
#     #     print('')
#     #     print('Pos:  ', template)
#     #     print('')
#     #     print('####################################################')
#
#     return((template, template_ID, template_pept))
