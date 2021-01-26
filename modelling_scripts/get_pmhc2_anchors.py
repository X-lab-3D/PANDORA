#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 16:00:32 2020

@author: rafaella
"""

import os
import sys
import numpy as np
from Bio.PDB.PDBParser import PDBParser
#from Bio.PDB.PDBParser import PDBParser


def get_anchors(pdbfile):
    
    #uses the function below  pocket_residues () that adapts the pocket residue to the PDB numbering
    pocket_M, pocket_N = pocket_residues(pdbfile)
    cutoff = 18
    
    #get all contacts between the chains with cutoff 18
    #contact file will be saved in a folder named contact_files 
    if not os.path.isdir('./contact_files'):
        os.mkdir ('./contact_files')
        
    if "contact-chainID_allAtoms" not in os.listdir('../modelling_scripts'):
        os.popen('g++ ../modelling_scripts/contact-chainID_allAtoms.cpp -o ../modelling_scripts/contact-chainID_allAtoms').read()
    os.popen('../modelling_scripts/contact-chainID_allAtoms %s %s > contact_files/all_contacts_%s.list' %(pdbfile, cutoff, pdbfile.split('/')[-1].split('.')[0])).read()
   
        
    with open( './contact_files/all_contacts_%s.list' %pdbfile.split('/')[-1].split('.')[0], 'r') as contacts:                             # Template contacts
          
            # vectors that will be used to count the frequency; each value represents the contact frequency between the 
            #previously defined pocket residues and the peptide residue corresponding to the vector index
            contact_1 = [0]*30
            contact_2 = [0]*30
            contact_3 = [0]*30
            contact_4 = [0]*30
            
            #read each line in the contacts file
            #check for distances between the MHCII chains and the peptide that are smaller than values specified in the pocket_M/N dictionaries
            for line in contacts:
                
               #the folllowing if statement checks the order of P, M in the contact file  
               chain_1 = line.split ("\t")[1]
               chain_2 = line.split("\t")[6]
               
               
               if chain_1 == 'P':
                    cp_aa_id = chain_1
                    cm_aa_id = chain_2
                    m_aa_id = line.split("\t")[7]  
                    p_aa_id = line.split("\t")[2]
                    atom_name = line.split ("\t")[3]
                    distance = line.split ("\t")[-1]   
               else: 
                   cm_aa_id = chain_1
                   cp_aa_id = chain_2
                   m_aa_id = line.split("\t")[2]
                   p_aa_id = line.split("\t")[7]
                   atom_name = line.split ("\t")[8]
                   distance = line.split ("\t")[-1]
                       
            
            #checks only the lines concerning the C alpha atom of the peptide residues/ C1 atom of modified residues
               if cp_aa_id == 'P':
                    if atom_name == ' CA ' or atom_name == ' C1 ':
                                
       
                        #Checks pocket in chain M
                        if cm_aa_id == 'M':
                            
                            n=0;
                            for aa in pocket_M.get('pocket_anch1'):
                                 
                               try:
                                       if int(m_aa_id) == aa:
                                           if (float(distance) <= pocket_M.get('pocket_1_distance')[n]):
                                               contact_1[int(p_aa_id)]+=1
                                          
                               except ValueError:
                                   continue
                               n+=1;
                               
                               
                            n=0;
                            for aa in pocket_M.get('pocket_anch3'):
    
                               try:
                                       if int(m_aa_id) == aa:
                                           if (float(distance) <= pocket_M.get('pocket_3_distance')[n]):
                                               contact_3[int(p_aa_id)]+=1
                                          
                               except ValueError:
                                   continue
                               n+=1;
                           
                                
                           
                        #Check pocket in chain N       
                        if cm_aa_id == 'N':
                            
                              
                            n=0;
                            for aa in pocket_N.get('pocket_anch2'):
                                 
                               try:
                                       if int(m_aa_id) == aa:
                                           if (float(distance) <= pocket_N.get('pocket_2_distance')[n]):
                                               contact_2[int(p_aa_id)]+=1
                                          
                               except ValueError:
                                   continue
                               n+=1;
                             
                             
                                
                            n=0;
                            for aa in pocket_N.get('pocket_anch3'):
                                 
                               try:
                                       if int(m_aa_id) == aa:
                                           if (float(distance) <= pocket_N.get('pocket_3_distance')[n]):
                                               contact_3[int(p_aa_id)]+=1
                                          
                               except ValueError:
                                   continue
                               n+=1;
                               
                            n=0;
                            for aa in pocket_N.get('pocket_anch4'):
                                 
                               try:
                                       if int(m_aa_id) == aa:
                                           if (float(distance)<= pocket_N.get('pocket_4_distance')[n]):
                                               contact_4[int(p_aa_id)]+=1                                            
                                          
                               except ValueError:
                                   continue
                               n+=1;
                            
 
                                              
    anchor_1 = np.argmax(contact_1) 
    anchor_2 = np.argmax(contact_2) 
    anchor_3 = np.argmax(contact_3) 
    anchor_4 = np.argmax(contact_4) 
    
    return (anchor_1, anchor_2, anchor_3, anchor_4)
    
  
    
###FUNCTION THAT DETERMINES POCKET
def pocket_residues (pdbfile):
    ### function takes as input pdbfile and outputs 2 dictionaries that define the peptide binding pocket in chains M, N of MHCII
    
    #normal pocket  numbering; pocket_anchn defines the residues of chain M/N that form the binding pocket of anchor n
    #pocket_n_distance represents the distances between the nth anchor residue and corresponding nth pocket residues from the first list
    pocket_M= {'pocket_anch1' : [28, 33, 31], 'pocket_1_distance' : [15, 10.5, 18], 
                'pocket_anch3' : [ 10], 'pocket_3_distance' : [13.5]}
    pocket_N = {'pocket_anch2' : [ 10, 11, 21, 22, 23], 'pocket_2_distance': [10.7, 10.5, 14.7, 10, 14],
               'pocket_anch3': [ 6, 7], 'pocket_3_distance': [13, 7.5],
               'pocket_anch4': [ 28, 27, 31, 33], 'pocket_4_distance': [14, 12.4, 15.7, 13]}
    
    #the pocket residues are defined differently when the numbering in the M or/and N chain starts from 1000
    pocket_M_renumbered = {'pocket_anch1' : [1028, 1033, 1031], 'pocket_1_distance' : [15, 10.5, 18], 
                           'pocket_anch3' : [ 1010], 'pocket_3_distance' : [13.5]}
    pocket_N_renumbered= {'pocket_anch2' : [ 1010, 1011, 1021, 1022, 1023], 'pocket_2_distance': [10.7, 10.5, 14.7, 10, 14],
                          'pocket_anch3': [1006, 1007], 'pocket_3_distance': [13, 7.5],
                          'pocket_anch4': [ 1028, 1027, 1031, 1033], 'pocket_4_distance': [14, 12.5, 15.7, 13]}
    
    
    ID = pdbfile.split('/')[-1].split('.')[0]
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure( ID, pdbfile)
    
    #for each chain (M/N), the first 5 residues are checked; if one of them is >1000, 
    #the pocket is defined using the pocket_M/N_renumbered dictionary
    
    for chain in structure.get_chains():
        if (chain.get_id() == 'M'):
            n = 0;
            for residue in chain.get_residues():
                n+=1;
                if n <= 7:
                    try:
                        if int (residue.id[1] > 1000):
                            pocket_M = pocket_M_renumbered
                    except ValueError:
                        continue
                else:
                    break;
                    
        if (chain.get_id() == 'N'):
            n = 0;
            for residue in chain.get_residues():
                n+=1;
                if n <= 7:
                    try:
                        if int (residue.id[1] > 1000):
                            pocket_N = pocket_N_renumbered
                    except ValueError:
                        continue
                else:
                    break;
                    
    return (pocket_M, pocket_N)
                
                    
           
    

  
