#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 00:50:19 2020

@author: rafaella
"""
import os
import sys
import numpy as np
from Bio.PDB.PDBParser import PDBParser



def small_molecule (pdbfile, chain_M, chain_P):
    #function takes as input a file (which can be either a pdb file or a contacts list), name of chain M and name of chain P
    #function outputs the name of the small molecule if there is one inside the binding pocket; otherwise, returns None
    
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
    #print (count_molecule)
    if np.max(count_molecule)>0:
        return molecule_in_pocket[np.argmax(count_molecule)]
    else:
        return None
               


#### checking the output:
def chain_peptide (pdbfile):
    
    #function that returns chain ID of peptide 
    
   ID = pdbfile.split('/')[-1].split('.')[0]
   parser = PDBParser(PERMISSIVE=1)
   structure = parser.get_structure( ID, pdbfile)
   chain_peptide = ''
   for chain in structure.get_chains():
        length = 0
        for residue in chain.get_residues():
            length +=1
        if length > 7 and length < 25:
           chain_peptide = chain.get_id()
           break;
   return chain_peptide

#%%
#ID = '1K5N'
#print ( type (small_molecule ('./data/PDBs_SM/%s.pdb' %ID , 'A', 'C') ))
#print ( small_molecule ('./data/PDBs_SM/%s_modified.pdb' %ID, 'A', 'C'))


#%%

#IDs = ['1XR8', '2GTW', '5FA4','5FA4_modified', '3UPR', '1QO3', '1K5N', '1K5N_modified', '1K5N_modified_2', '1JGD', #'1JGD_modified', '2ESV','2ESV_modified', '3DX6','3RWG', '4NNX', '5FDW', '6O4Z']
#with open("small_molec.txt", 'w') as output:
#    for ID in IDs:
#        pdbfile ='./data/PDBs_SM/%s.pdb' %ID
#        contact = './all_contacts_%s.list' %ID
#        chain_pep = chain_peptide(pdbfile)
#        print (pdbfile.split('/')[-1].split('.')[0], small_molecule (pdbfile, 'A', '%s' %chain_pep), file = output)