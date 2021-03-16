from PANDORA.Contacts import Contacts
import numpy as np


### ### Function that finds the anchors
def get_anchors_pMHCI(pdb):
    ''' Find the peptide anchor residues of a pMHCI structure using the Contacts class

    Args:
        pdb: Bio.PDB object
    Returns: tuple of anchor residues (int)

    '''

    cutoff = 11.2
    ###Define pocket residues
    pocket = { 'pocket_anch1' : [24, 25, 35, 36],
              'pocket_anch2' : [81, 1005, 1027, 1028, 1029, 1030, 1033]}
    pocket_wrong_numbering = { 'pocket_anch1' : [1024, 1025, 1035, 1036],
                    'pocket_anch2' : [1081, 2005, 2027, 2028, 2029, 2030, 2033]}

    # Template contacts
    contacts = Contacts.Contacts(pdb, cutoff=cutoff).chain_contacts

    contact_1 = [0]*25
    contact_2 = [0]*25

    #for cases where numbering starts from 1000, we use the updated pocket residues
    if int(contacts[0][6]) >1000:
        pocket = pocket_wrong_numbering

    for line in contacts:
        #looks at one residue at a time
        m_aa_id = line[6]
        p_aa_id = line[2]
        atom_name = line[3]

        if atom_name == 'CA':
            for aa in pocket.get('pocket_anch1'):
                try:
                    if int(m_aa_id) == aa:
                        contact_1[int(p_aa_id)]+=1
                except ValueError:
                    continue

            for aa in pocket.get('pocket_anch2'):
                try:
                    if int(m_aa_id) == aa:
                        contact_2[int(p_aa_id)]+=1
                except ValueError:
                    continue

    anchor_1 = np.argmax(contact_1)
    anchor_2 = np.argmax(contact_2)

    return [anchor_1, anchor_2]


def pMHCII_anchors(pdb):
    ''' Finds peptide anchor residues of p:MHCII complexes

    :param pdb: Bio.PDB object (or path to pdb file (string))
    :return: (list of ints) Peptide anchor residue positions
    '''

    # Rafaella's residues
    pocket_M = {'anch1':[28, 33, 31], 'anch3':[10]}
    pocket_N = {'anch2':[10, 11, 21, 22, 23], 'anch3':[6,7], 'anch4':[ 28, 27, 31, 33]}

    # Derek's slightly modified version
    # pocket_M = {'anch1':range(28,36), 'anch3':[10]}
    # pocket_N = {'anch2':[10, 11, 21, 22, 23], 'anch3':[6,7], 'anch4':range(33,37)}


    # Calculate the contacts with a cutoff of 18. This cutoff because at lower cutoffs, the orientations of side chains
    # of different structures, result in highly different contacts.
    cont = Contacts.Contacts(pdb, cutoff=18).chain_contacts

    # Make sure the first chain is P, second chain is M or N, not the other way around
    cont2 = []
    for i in cont:
        if i[1] == 'P' or i[5] == 'P':
            cont2.append(i)
            cont2.append((i[4], i[5], i[6], i[7], i[0], i[1], i[2], i[3], i[8]))
    cont = list(dict.fromkeys(cont2))
    del(cont2)

    # Find the peptide residues that reside in the M and N pockets
    anch1, anch2, anch3, anch4 = [], [], [], []
    for i in cont:
        if i[1] == 'P' and i[5] == 'M':
            if i[3] != 'CA' or i[3] != 'CB' or i[3] != 'CG':
                if i[6] in pocket_M['anch1']:
                    anch1.append(i)
                # if i[6] in range(7,9): # 8,15
                #     anch2.append(i)
                if i[6] in pocket_M['anch3']:
                    anch3.append(i)

        if i[1] == 'P' and i[5] == 'N':
            if i[3] != 'CA' or i[3] != 'CB' or i[3] != 'CG':
                if i[6] in pocket_N['anch2']: #10 11 21, 22, 23          mijn: 13, 14, 15, 16, 17
                    anch2.append(i)
                if i[6] in pocket_N['anch3']: #6 7
                    anch3.append(i)
                if i[6] in pocket_N['anch4']:  # [28, 27, 31, 33]   # Mijn: range 33 37
                    anch4.append(i)

    # Find the peptide residue with the minimal mean distance to the specified pocket residues.
    try:
        anch1_list = []
        for i in set([i[2] for i in anch1]):
            l = [x[-1] for x in anch1 if x[2] == i]
            anch1_list.append((i, sum(l)/len(l) ))
        anchor_1 = min((x[1], x[0]) for x in anch1_list)[1]
    except:
        anchor_1 = 0
    try:
        anch2_list = []
        for i in set([i[2] for i in anch2]):
            l = [x[-1] for x in anch2 if x[2] == i]
            anch2_list.append((i, sum(l)/len(l) ))
        anchor_2 = min((x[1], x[0]) for x in anch2_list)[1]
    except:
        anchor_2 = 0
    try:
        anch3_list = []
        for i in set([i[2] for i in anch3]):
            l = [x[-1] for x in anch3 if x[2] == i]
            anch3_list.append((i, sum(l)/len(l) ))
        anchor_3 = min((x[1], x[0]) for x in anch3_list)[1]
    except:
        anchor_3 = 0
    try:
        anch4_list = []
        for i in set([i[2] for i in anch4]):
            l = [x[-1] for x in anch4 if x[2] == i]
            anch4_list.append((i, sum(l)/len(l) ))
        anchor_4 = min((x[1], x[0]) for x in anch4_list)[1]
    except:
        anchor_4 = 0

    # todo properly test this. I changed Rafaella's pocket positions slightly in some cases. Needs to be optimized.
    # Anchor 2 and 3 can be difficult. If anchor 1 and 4 are spaced in the canonical way, and anchor 2 and 3 are not,
    # then just jump the gun and use the canonical spacing for anchor 2 and 3.
    if anchor_4 - anchor_1 == 8:
        if anchor_2 - anchor_3 != -2 or anchor_3 - anchor_2 != 2:
            anchor_2 = anchor_1 + 3
            anchor_3 = anchor_1 + 5

    # Return the ordered set of anchors. Filter out 0's and >25 numbers (don't know why these can be in here though)
    return sorted(set([i for i in [anchor_1, anchor_2, anchor_3, anchor_4] if i != 0 or i > 25]))