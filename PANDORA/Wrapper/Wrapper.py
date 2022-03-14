#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Wrapper.run_model import run_model
import csv
from joblib import Parallel, delayed

class Wrapper():
    def __init__(self):

        """Pandora wrapper object.

        Args:
            None.

        Returns:
            None.

        """
        self.data_file = ''
        self.db = None
        self.targets = {}
        self.jobs = {}


    def __get_targets_from_file(self, data_file, delimiter='\t', header=True,
                               IDs_col=None, peptides_col=0,
                               allele_col=1, anchors_col=None,
                               M_chain_col=None, N_chain_col=None,
                               start_row=None, end_row=None):
        """
        Extracts peptide sequences, alleles and anchors (if specified)
            from the target file.
            Default input should be a .tsv file without any header with
            the following structure: peptides_sequence_col \t alleles_name_col

        Args:
            data_file (str): Path to the input tsv/csv file containing targets
                information.
            delimiter (str, optional): data_file delimiter. Do not use
                semicolons (';') as separators. Defaults to '\t'.
            header (bool, optional): If True, assumes the data_file has a
                header line and skips it. If your file has no header line,
                set it as False. Defaults to True.
            IDs_col (int or None, optional): Column of data_file containing
                the targets IDs. If None, will automatically assign an ID
                according to the row number. Defaults to None.
            peptides_col (int, optional): Column of data_file containing
                the targets peptides. Defaults to 0.
            allele_col (int, optional): Column of data_file containing
                the targets alleles. Umbiguous allele cases (where the allele
                might have multiple names) should be separated by a
                semicolon (';'). Defaults to 1.
            anchors_col (None or int, optional): Column of data_file containing
                the targets anchors. Anchors should be two numbers separated
                by a semicolon (';'). Defaults to 2.
            M_chain_col (None or int, optional): Column of data_file containing
                the targets M chain sequences.
            N_chain_col (None or int, optional): Column of data_file containing
                the targets N chain sequences (only for MHCII).
            start_row (None or int, optional): Starting row of data_file, to use when
                splitting the data_file into multiple batches. This allows to
                specify from which row the samples for this job start.
            end_row (None or int, optional): Ending row of data_file, to use when
                splitting the data_file into multiple batches. This allows to
                specify at which row the samples for this job end.

        Returns:
            None.

        """

        targets = {}
        with open(data_file, 'r') as infile:
            spamreader = csv.reader(infile, delimiter=delimiter)
            if header == True:
                next(spamreader)
            for i, row in enumerate(spamreader):
                if start_row != None and i < start_row:
                    pass
                elif end_row != None and i >= end_row:
                    break
                else:
                    ## Assign target ID
                    if IDs_col != None:
                        target_id = row[IDs_col]
                    else:
                        target_id = 'Target_%i' %(i+1)

                    ## Assign peptide sequence
                    peptide_seq = row[peptides_col]

                    ## Assign allele name
                    allele = row[allele_col].split(';')

                    ## Make target entry
                    targets[target_id] = {'peptide_sequence' : peptide_seq,
                                              'allele' : allele}

                    ## Assign optional arguments. Be sure the empty values correspond
                    ## to the default values in PMHC.Target.__init__()
                    ## Assign anchors
                    if anchors_col:
                        anchors = tuple([int(x) for x in row[anchors_col].split(';')])
                        targets[target_id]['anchors'] = anchors
                    else:
                        targets[target_id]['anchors'] = []

                    ## Assign M chain sequence
                    if M_chain_col:
                        M_chain_seq = row[M_chain_col]
                        targets[target_id]['M_chain_seq'] = M_chain_seq
                    else:
                        targets[target_id]['M_chain_seq'] = ''

                    ## Assign N chain sequence
                    if N_chain_col:
                        N_chain_seq = row[N_chain_col]
                        targets[target_id]['N_chain_seq'] = N_chain_seq
                    else:
                        targets[target_id]['N_chain_seq'] = ''

        self.targets = targets

    def create_targets(self, data_file, database, MHC_class, delimiter = '\t',
                       header=True, IDs_col=None, peptides_col=0, allele_col=1,
                       anchors_col=None, M_chain_col=None, N_chain_col=None,
                       benchmark=False, verbose=False, start_row=None,
                       end_row=None, use_netmhcpan=False):
        """


        Args:
            data_file (str): Path to the input tsv/csv file containing targets
                information.
            database (PANDORA.Database.Database): Database object.
            MHC_class (str): MHC class of the targets, as 'I' or 'II'.
            delimiter (str, optional): data_file delimiter. Do not use
                semicolons (';') as separators. Defaults to '\t'.
            header (bool, optional): If True, assumes the data_file has a
                header line and skips it. If your file has no header line,
                set it as False. Defaults to True.
            IDs_col (int or None, optional): Column of data_file containing
                the targets IDs. If None, will automatically assign an ID
                according to the row number. Defaults to None.
            peptides_col (int, optional): Column of data_file containing
                the targets peptides. Defaults to 0.
            allele_col (int, optional): Column of data_file containing
                the targets alleles. Umbiguous allele cases (where the allele
                might have multiple names) should be separated by a
                semicolon (';'). Defaults to 1.
            anchors_col (int, optional): Column of data_file containing
                the targets anchors. Anchors should be two numbers separated
                by a semicolon (';'). Defaults to 2.
            M_chain_col (None or int, optional): Column of data_file containing
                the targets M chain sequences.
            N_chain_col (None or int, optional): Column of data_file containing
                the targets N chain sequences (only for MHCII).
            benchmark (bool, optional): Set True only for benchmarking purpose,
                if target structures are available. Defaults to False.
            start_row (None or int): Starting row of data_file, to use when
                splitting the data_file into multiple batches. This allows to
                specify from which row the samples for this job start.
            end_row (None or int): Ending row of data_file, to use when
                splitting the data_file into multiple batches. This allows to
                specify at which row the samples for this job end.
            use_netmhcpan (bool, optional): If True, uses local installation
                of netMHCPan to predict anchor positions for each target.

        Returns:
            None.

        """
        self.data_file = data_file
        self.db = database

        ## Extract targets from data_file
        self.__get_targets_from_file(data_file, delimiter=delimiter,
                                     header=header, IDs_col=IDs_col,
                                     peptides_col=peptides_col, allele_col=allele_col,
                                     anchors_col=anchors_col, M_chain_col=M_chain_col,
                                     N_chain_col=N_chain_col, start_row=start_row,
                                     end_row=end_row)

        ## Create target objects
        jobs = {}
        for target_id in self.targets:
            #try:
                if verbose:
                    print('Target ID: ', target_id)
                    print('Target MHC_class: ', MHC_class)
                    print('Target allele: ', self.targets[target_id]['allele'])
                    print('Target peptide: ', self.targets[target_id]['peptide_sequence'])
                    print('Target M chain seq: ', self.targets[target_id]['M_chain_seq'])
                    if N_chain_col:
                        print('Target N chain seq: ', self.targets[target_id]['N_chain_seq'])
                if verbose:
                    print('Target Anchors: ', self.targets[target_id]['anchors'])
                #try:
                tar = PMHC.Target(target_id, allele_type=self.targets[target_id]['allele'],
                                  peptide=self.targets[target_id]['peptide_sequence'] ,
                                  MHC_class=MHC_class, anchors=self.targets[target_id]['anchors'],
                                  M_chain_seq=self.targets[target_id]['M_chain_seq'],
                                  N_chain_seq=self.targets[target_id]['N_chain_seq'],
                                  use_netmhcpan=use_netmhcpan)
                #except Exception as err:
                #    print('Skipping Target %s at Target object generation step for the following reason:' %target_id)
                #    print(("Exception: {0}".format(err)))
                try:
                    mod = Pandora.Pandora(tar, self.db)
                except Exception as err:
                    print('Skipping Target %s at Pandora object generation step for the following reason:' %target_id)
                    print(("Exception: {0}".format(err)))
                try:
                    mod.find_template(benchmark=benchmark)
                    jobs[target_id] = [tar, mod.template]
                except Exception as err:
                    print('Skipping Target %s at template selection step for the following reason:' %target_id)
                    print(("Exception: {0}".format(err)))
            #except Exception as err:
            #    print('An unidentified problem occurred with Target %s. Please check your target info' %target_id)
            #    print(("Exception: {0}".format(err)))
        self.jobs = jobs

    def run_pandora(self, num_cores=1, n_loop_models=20, n_jobs=None,
                    benchmark=False, output_dir=False, pickle_out=False):
        """Runs Pandora in parallel jobs.


        Args:
            num_cores (int, optional): Number of parallel PANDORA jobs.
                Each one will be sent to a different core. Defaults to 1.
            n_loop_models (int, optional): Number of  loop models.
                Defaults to 20.
            n_jobs (int, optional): Number of parallel MODELLER loop jobs.
                Do not increase further than n_loop_models. Defaults to None.
            benchmark (bool, optional): Set True only for benchmarking purpose,
                if target structures are available. Defaults to False.
            output_dir (str, optional): Output directory path.
                Defaults to False.
            pickle_out (bool, optional): If True, outputs a pickle file
                containing every model object. Defaults to False.

        Returns:
            None.

        """

        for job in self.jobs:
            if output_dir:
                self.jobs[job].extend([n_loop_models, n_jobs, benchmark, pickle_out, output_dir])
            else:
                self.jobs[job].extend([n_loop_models, n_jobs, benchmark, pickle_out])
        Parallel(n_jobs = num_cores, verbose = 1)(delayed(run_model)(job) for job in list(self.jobs.values()))
