#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 14:05:06 2020

@author: Rafaella Buzatu
"""

import os
#import sys
import numpy as np
## Inputs: pdb file  Outputs: anchor positions

### ### Function that finds the anchors
def get_anchors(pdbfile, rm_outfile=False):

    cutoff =11.2
    ###Define pocket residues
    pocket = { 'pocket_anch1' : [ 24, 25, 35, 36],
              'pocket_anch2' : [81, 1005, 1027, 1028, 1029, 1030, 1033]}
    pocket_wrong_numbering = { 'pocket_anch1' : [1024, 1025, 1035, 1036],
                    'pocket_anch2' : [1081, 2005, 2027, 2028, 2029, 2030, 2033]}


    ### Calculating all Atom contacts
    os.popen('../../../modelling_scripts/contact-chainID_allAtoms %s %s > ./all_contacts_for_anchors_%s.list' %(pdbfile, cutoff, pdbfile.split('/')[-1].split('.')[0])).read()
    
    with open( './all_contacts_for_anchors_%s.list' %pdbfile.split('/')[-1].split('.')[0], 'r') as contacts:                             # Template contacts
       # with open('./anchors_files/peptide_anchors_%s.list' %pdbfile.split('/')[-1].split('.')[0], 'w') as output:
    
            contact_1 = [0]*13
            contact_2 = [0]*13
          
            
            #for cases where numbering starts from 1000, we use the updated pocket residues
            line1=contacts.readline()
            m_aa_id_check = line1.split("\t")[2]   
            if int(m_aa_id_check) >1000:
                        pocket = pocket_wrong_numbering 
                
            for line in contacts:
               #looks at one residue at a time
                m_aa_id = line.split("\t")[2]    
                p_aa_id = line.split("\t")[7]
                atom_name = line.split ("\t")[8]
        
                
                if atom_name == ' CA ':    
                    
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
           
    if rm_outfile == True:
        os.popen('rm ./all_contacts_for_anchors_%s.list' %pdbfile.split('/')[-1].split('.')[0])
          
    return (anchor_1, anchor_2)
    #return (contact_1, contact_2)
