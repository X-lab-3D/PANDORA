#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 00:50:19 2020

@author: rafaella
"""
import os
import sys
import numpy as np
import pickle
from Bio.PDB.PDBParser import PDBParser



def small_molecule (pdbfile, chain_M, chain_P):
    #function takes as input a file (which can be either a pdb file or a contacts list), name of chain M and name of chain P
    #function outputs the name of the small molecule if there is one inside the binding pocket; otherwise, returns None
    
    print(pdbfile)
    cutoff =6
    pocket = [ 7, 9, 26, 28, 1007, 1009, 1024, 1026]


    ### Calculating all Atom contacts; 
    
    #if the input is a pdb file, the contacts file is calculated
    if pdbfile.endswith('.pdb'):
        contactfile = 'all_contacts_%s.list' %pdbfile.split('/')[-1].split('.')[0]
        os.popen('./modelling_scripts/contact-chainID_allAtoms %s %s > %s ' %( pdbfile, cutoff, contactfile)).read()
    
    #if the contacts file list is used as input, it is used further in the script
    elif pdbfile.endswith('.list'):
        contactfile = pdbfile

    else:
        return None
        
    with open(contactfile) as contacts:
                              
            molecule_in_pocket = []  #list where the name of small molecule(s) is/are appended when found in contact list
            count_molecule = [ 0, 0, 0, 0, 0]  #counter for the frequency of contacts between the small moelcule of corresponding index and the pocket
         
            
            #for cases where numbering starts from 1000, we use the updated pocket residues
            line1=contacts.readline()
            m_aa_id_check = line1.split("\t")[2]   
            if int(m_aa_id_check) >1000:
                for i in pocket:
                    i = i+1000
              
                
            for line in contacts:
                cm_aa_id = line.split("\t")[1]
                m_aa_id = line.split("\t")[2]    
                cp_aa_id = line.split("\t")[6]
                molecule = line.split("\t")[5]
        
                #check for contacts that are not between the MHC chain and peptide or water molecules
                if cm_aa_id == chain_M:
                    if cp_aa_id != chain_P  and cp_aa_id != 'B' and molecule != 'HOH':
                        
                        #once a suitable contact is found, molecule name is appended in moelcule_in_pocket if not previously there
                        #and the frequency count increases
                        for aa in pocket:
                           try:
                                   if int(m_aa_id) == aa:
                                       
                                       if molecule in molecule_in_pocket:
                                           count_molecule [molecule_in_pocket.index(molecule)] += 1
                                       else:   
                                            molecule_in_pocket.append(molecule);
                                            count_molecule [molecule_in_pocket.index(molecule)] = 1
                                                          
                           except ValueError:
                               continue
       
    #print (np.argmax(count_molecule))
    print (count_molecule)
    os.system('rm %s' %contactfile)
    if np.max(count_molecule)>0:
        return molecule_in_pocket[np.argmax(count_molecule)]
    else:
        return None
               


#### checking the output:
outputs = {}
for f in os.listdir('./data/PDBs/pMHCI/'):
    out = small_molecule('./data/PDBs/pMHCI/' + f, 'M', 'P')
    if out != None:
        outputs[f.split('_')[0]] = out

outfile = open('./data/csv_pkl_files/small_molecule_templates.pkl', 'wb')
pickle.dump(outputs, outfile)
outfile.close()
