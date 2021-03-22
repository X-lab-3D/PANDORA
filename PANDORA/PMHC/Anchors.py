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

    Args:
        pdb: Bio.PDB object (or path to pdb file (string))

    Returns: (list of ints) Peptide anchor residue positions

    '''

    def closest_pept_res(dist_list):
        '''Find the closest peptide residue based on te minimal mean value between peptide residues and MHC residues'''
        anch_list = []
        for i in set([i[2] for i in dist_list]):
            l = [x[-1] for x in dist_list if x[2] == i]
            anch_list.append((i, sum(l) / len(l))) #add mean distance between pept residue and MHC residue
        return min(((x[1], x[0]) for x in anch_list), default=(0, 0))[1] #get the pept residue with the min mean dist

    # Rafaella's residues
    # pocket_M = {'anch1':[28, 33, 31], 'anch3':[10]}
    # pocket_N = {'anch2':[10, 11, 21, 22, 23], 'anch3':[6,7], 'anch4':[ 28, 27, 31, 33]}

    # Derek's modified version. Empirically tested to have the highest prediction accuracy. See below. (only best shown)
    pocket_M = {'anch1':[28,29], 'anch4':[82,81,80], 'anch3':[113,114]}
    pocket_N = {'anch2':[24,25]}

    # Accuracies peptide residue anchor prediction using MHCII residue x without extra rules. Tested on all IMGT MHCII's
    # anchor 1: M28: 0.971, M29: 0.961
    # anchor 2: N24: 0.962, N25: 0.923
    # anchor 3: M113: 0.914, M114: 0.904
    # anchor 4: M80: 0.804, M81: 0.843, M81: 0.843

    # With extra rules (e.g. anchor 4 can only be n spaces from anchor 2 etc...)
    # Getting all anchors right --> Accuracy: 97.0588
    # Anchor 1: 98.0392
    # Anchor 2: 98.0392
    # Anchor 3: 99.0196
    # Anchor 4: 99.0196

    # Calculate the contacts with a cutoff of 18. This cutoff because at lower cutoffs, the orientations of side chains
    # of different structures, result in highly different contacts.
    cont = Contacts.Contacts(pdb, pept_contacts=True, M_only=sum(pocket_M.values(), []),
                             N_only=sum(pocket_N.values(), []), cutoff=25).chain_contacts

    # Make sure the first chain is P, second chain is M or N, not the other way around
    cont = cont + [(i[4], i[5], i[6], i[7], i[0], i[1], i[2], i[3], i[8]) for i in cont if i[1] == 'P' or i[5] == 'P']
    cont = list(dict.fromkeys(cont))

    # Find the peptide residues that reside in the M and N pockets
    anch1, anch2, anch3, anch4 = [], [], [], []
    for i in cont:
        if i[1] == 'P' and i[5] == 'M': # Use residues from the M chain to find anchor 1, 3 and 4
            # Check the Alpha, Beta and Gamma carbon molecules of the peptide residues.
            if i[3] == 'CA' or i[3] == 'CB' or i[3] == 'CG':
                if i[6] in pocket_M['anch1']: #If connected (<25 ang) to specified anch 1 MHC residues, add to list
                    anch1.append(i)
                if i[6] in pocket_M['anch3']: #If connected (<25 ang) to specified anch 3 MHC residues, add to list
                    anch3.append(i)
                if i[6] in pocket_M['anch4']: #If connected (<25 ang) to specified anch 4 MHC residues, add to list
                    anch4.append(i)


        if i[1] == 'P' and i[5] == 'N': # Use chain N for the second anchor
            # Check the Alpha, Beta and Gamma carbon molecules of the peptide residues
            if i[3] == 'CA' or i[3] == 'CB' or i[3] == 'CG':
                if i[6] in pocket_N['anch2']: #If connected (<25 ang) to specified anch 2 MHC residues, add to list
                    anch2.append(i)

    # Now find the peptide residues that are closest to the specified MHCII residues to get the anchors

    # First check if the peptide is inverted: C to N term instead of N to C term (this is very rare)
    anchor_1 = closest_pept_res(anch1)
    anchor_4 = closest_pept_res(anch4)
    # Is anchor 1 smaller than anchor 4? Yes --> peptide is inverted.
    if anchor_4 < anchor_1:
        anchor_2 = closest_pept_res(anch2)
        anchor_3 = closest_pept_res(anch3)

    # If its not, use some extra rules to increase prediction accuracy.
    else:
        # First find anchor 2, this anchor is found with the highest accuracy
        anchor_2 = closest_pept_res(anch2)

        # Make sure anchor 1 if smaller than anchor 2, but bigger than anchor 2 - 5
        # This makes sure residues from long flanking regions that are flipped over the MHC structure aren't selected
        anch1 = [i for i in anch1 if i[2] < anchor_2 - 1 and i[2] > anchor_2 - 5]
        anchor_1 = closest_pept_res(anch1)

        # Make sure anchor 3 is close to anchor 2 as well (but not too far, to prevent selecting anchor 4)
        anch3 = [i for i in anch3 if i[2] > anchor_2 + 1 and i[2] < anchor_2 + 4]
        anchor_3 = closest_pept_res(anch3)

        # Same for anchor 4. Make sure its not too far away. Otherwise distant flanking residues might get selected.
        anch4 = [i for i in anch4 if i[2] > anchor_2 + 2 and i[2] < anchor_2 + 6]
        anchor_4 = closest_pept_res(anch4)

    # Put the anchors in the right order if they aren't already
    anchors = sorted([i for i in [anchor_1, anchor_2, anchor_3, anchor_4]])

    # Anchor 3 can be difficult in extreme cases. If anchor 1 and 4 are spaced in the canonical way,
    # and anchor 2 and 3 are not, then just jump the gun and use the canonical spacing for anchor 2 and 3.
    if anchors[-1] - anchors[0] == 8:
        if anchors[1] - anchors[2] != -2 or anchors[2] - anchors[1] != 2:
            anchors[1] = anchors[0] + 3
            anchors[2] = anchors[0] + 5

    # Return the ordered set of anchors. Filter out 0's and >25 numbers if they happen to get through.
    return sorted(set([i for i in anchors if i != 0 or i > 25]))
