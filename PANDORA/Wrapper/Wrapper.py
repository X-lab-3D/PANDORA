#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 17:46:39 2021

@author: Dario Marzella
"""
import PANDORA
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database
from PANDORA.Wrapper.run_model import run_model
from pathos.multiprocessing import ProcessingPool as Pool
from pathos.multiprocessing import freeze_support
import time
import csv
import os

#test joblib
from joblib import Parallel, delayed
from multiprocessing import Manager
from joblib.externals.loky import set_loky_pickler
set_loky_pickler("dill")


class Wrapper():
    def __init__(self):

        """

        Args:
            data_file (str): path to the csv/tsv file
            db (Database): PANDORA.Database object

        Returns:
            None.

        """
        self.data_file = ''
        self.db = None
        self.targets = {}
        self.jobs = {}
        
        
    def __get_targets_from_file(self, data_file, delimiter='\t', header=True, IDs_col=None, 
                                peptides_col=0, allele_col=1, anchors_col=None):
        """Extracts peptide sequences, alleles and anchors (if specified) from the target file.
           Default input should be a .tsv file without any header with the following structure:
               peptides_sequence_col \t alleles_name_col
        

        Args:
            data_file (str): path to .tsv or .csv file containing target sequences.
            delimiter (TYPE, optional): DESCRIPTION. Defaults to '\t'.
            header (TYPE, optional): DESCRIPTION. Defaults to True.
            IDs_col (TYPE, optional): DESCRIPTION. Defaults to None.
            peptides_col (TYPE, optional): DESCRIPTION. Defaults to 0.
            allele_col (TYPE, optional): DESCRIPTION. Defaults to 1.
            anchors_col (TYPE, optional): DESCRIPTION. Defaults to None.

        Returns:
            None.

        """

        targets = {}
        with open(data_file, 'r') as infile:
            spamreader = csv.reader(infile, delimiter=delimiter)
            if header == True:
                next(spamreader)
            for i, row in enumerate(spamreader):
                ## Assign target ID
                if IDs_col != None:
                    target_id = row[IDs_col]
                else:
                    target_id = 'Target_%i' %(i+1)
                    
                ## Assign peptide sequence
                peptide_seq = row[peptides_col]
                
                ## Assign allele name
                allele = row[allele_col].split(';')
                
                ## Assign anchors
                if anchors_col:
                    anchors = tuple([int(x) for x in row[anchors_col].split(',')])
                #if 'HLA' in allele:
                #    star_allele = (allele[0:5]+'*'+allele[5:])
                #    targets.append((target_id, seq, star_allele))
                #else:
                
                    targets[target_id] = {'peptide_sequence' : peptide_seq, 'allele' : allele,
                                          'anchors' : anchors}
                else:
                    targets[target_id] = {'peptide_sequence' : peptide_seq, 'allele' : allele}
        self.targets = targets

    def create_targets(self, data_file, db, MHC_class, delimiter = '\t', header=True, 
                       IDs_col=None, peptides_col=0, allele_col=1, anchors_col=None, 
                       benchmark=False, num_models=20, verbose=False):
        """
        

        Args:
            data_file (TYPE): DESCRIPTION.
            db (TYPE): DESCRIPTION.
            MHC_class (TYPE): DESCRIPTION.
            delimiter (TYPE, optional): DESCRIPTION. Defaults to '\t'.
            header (TYPE, optional): DESCRIPTION. Defaults to True.
            IDs_col (TYPE, optional): DESCRIPTION. Defaults to None.
            peptides_col (TYPE, optional): DESCRIPTION. Defaults to 0.
            allele_col (TYPE, optional): DESCRIPTION. Defaults to 1.
            anchors_col (TYPE, optional): DESCRIPTION. Defaults to 2.

        Returns:
            None.

        """
        self.data_file = data_file
        self.db = db
        
        ## Extract targets from data_file
        self.__get_targets_from_file(data_file, delimiter=delimiter, header=header, IDs_col=IDs_col, 
                                     peptides_col=peptides_col, allele_col=allele_col, anchors_col=anchors_col)
        
        ## Create target objects
        jobs = {}
        for target_id in self.targets:
            try:
                if verbose:
                    print('Target ID: ', target_id)
                    print('Target MHC_class: ', MHC_class)
                    print('Target allele: ', self.targets[target_id]['allele'])
                    print('Target peptide: ', self.targets[target_id]['peptide_sequence'])
                try:
                    if verbose:
                        print('Target Anchors: ', self.targets[target_id]['anchors'])
                    tar = PMHC.Target(target_id, allele_type=self.targets[target_id]['allele'],
                                      peptide=self.targets[target_id]['peptide_sequence'] ,
                                      MHC_class=MHC_class, anchors=self.targets[target_id]['anchors'])
                except KeyError:
                    tar = PMHC.Target(target_id, allele_type=self.targets[target_id]['allele'],
                                      peptide=self.targets[target_id]['peptide_sequence'],
                                      MHC_class=MHC_class)
                mod = Pandora.Pandora(tar, db)
                try:
                    mod.find_template(benchmark=benchmark)
                    jobs[target_id] = [tar, mod.template]
                except Exception:
                    print('Skipping Target %s' %target_id)
            except: ### TODO: test and specify exception for this except
                print('An unidentified problem occurred with Target %s. Please check your target info' %target_id)
        self.jobs = jobs
        
    def __run_multiprocessing(self, func, num_cores):
        """
        

        Args:
            func (function): DESCRIPTION.
            num_cores (int): DESCRIPTION.

        Returns:
            TYPE: DESCRIPTION.

        """
        with Pool(processes=num_cores) as pool:
            return pool.map(func, list(self.jobs.values()))

    def run_pandora(self, num_cores=1, num_models=20, n_jobs=None, benchmark=False):
        """
        

        Args:
            n_cores (int, optional): DESCRIPTION. Defaults to 1.
            benchmark (TYPE, optional): DESCRIPTION. Defaults to False.

        Returns:
            None.

        """
        
        freeze_support()
        for job in self.jobs:
            self.jobs[job].extend([num_models, n_jobs, benchmark])
        self.__run_multiprocessing(run_model, num_cores)
    
    def run_pandora_joblib(self, num_cores=1, num_models=20, n_jobs=None, benchmark=False):
        """
        

        Args:
            n_cores (TYPE, optional): DESCRIPTION. Defaults to 1.
            benchmark (TYPE, optional): DESCRIPTION. Defaults to False.

        Returns:
            None.

        """

        for job in self.jobs:
            self.jobs[job].extend([num_models, n_jobs, benchmark])
        Parallel(n_jobs = num_cores, verbose = 1, backend='loky')(delayed(run_model)(job) for job in list(self.jobs.values()))
'''
## To test

import dill
with open('PANDORA_files/data/csv_pkl_files/database.pkl', 'rb') as inpkl:
    db = dill.load(inpkl)
wrap = Wrapper()
wrap.create_targets('PANDORA_files/data/csv_pkl_files/test_datafile.tsv', db, 
                    MHC_class='I', header=False, IDs_col=0, peptides_col=1, 
                    allele_col=2, benchmark=True)
wrap.run_pandora(num_cores=6, num_models=20, benchmark=True)
MHC_class='I'
for target in wrap.targets:
    tar = PMHC.Target(target, allele_type=wrap.targets[target]['allele'], peptide=wrap.targets[target]['peptide_sequence'], MHC_class=MHC_class)
    break
tar.allele_type
mod = Pandora.Pandora(tar, db)
mod.find_template()
'''
    

