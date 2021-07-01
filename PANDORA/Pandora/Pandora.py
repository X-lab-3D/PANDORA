
import PANDORA
from PANDORA.Pandora import Align
from PANDORA.Pandora import Modelling_functions
import time
import os
from Bio.PDB import PDBParser



class Pandora:

    def __init__(self, target, database = None, template = None, output_dir = PANDORA.PANDORA_data + '/outputs'):
        self.target = target
        self.template = template
        self.database = database
        self.output_dir = output_dir

        if database == None and template == None:
            raise Exception('Provide a Database object so Pandora can find the best suitable template structure for '
                            'modelling. Alternatively, you can specify a user defined Template object.')

    def find_template(self, seq_based_templ_selection=False, benchmark=False, verbose=True):
        ''' Find the best template structure given a Target object

        Args:
            seq_based_templ_selection: (bool) Use template selection based on template sequences instead of allele.
            verbose: (bool) Print information

        '''

        if verbose:
            print('\tTarget MHC Class: %s' % self.target.MHC_class)
            print('\tTarget Allele:  %s' % self.target.allele_type)
            print('\tTarget Peptide: %s' % self.target.peptide)
            print('\tTarget Anchors: %s,%s\n' % (self.target.anchors[0],self.target.anchors[1]))

        if self.template == None: # Only find the best template if the user didn't specify one
            if verbose and self.target.M_chain_seq != '' and seq_based_templ_selection:
                print('\tUsing sequence based template selection')
            elif verbose:
                print('\tUsing allele type based template selection')
            # Find the best template. If the target already exists in the database, 
            # also consider the initial loop model as a model
            self.template, self.keep_IL = Modelling_functions.find_template(self.target, self.database,
                                                              seq_based_templ_selection=seq_based_templ_selection,
                                                              benchmark=benchmark)
            self.target.templates = [self.template.id]
            if verbose:
                print('\tSelected template structure: %s' %self.template.id)

        else:
            if verbose:
                print('\tUser defined template structure: %s' %self.template.id)
            # Check if the target structure and template structure are the same.
            self.keep_IL = Modelling_functions.check_target_template(self.target, self.template)

        if verbose:
            print('\tTemplate Allele:  %s' % self.template.allele_type)
            print('\tTemplate Peptide: %s' % self.template.peptide)
            print('\tTemplate Anchors: %s,%s\n' % (self.target.anchors[0],self.target.anchors[1]))

    def prep_output_dir(self):
        ''' Create an output directory and move the template pdb there
        '''
        # create an output directory
        try:
            self.output_dir = '%s/%s_%s' %(self.output_dir, self.target.id , self.template.id)
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
        except:
            pass
        os.system('cp %s %s/%s.pdb' %(self.template.pdb_path, self.output_dir, self.template.id))

    def align(self, verbose = True):
        ''' Create the alignment file for modeller

        Args:
            verbose:  (bool) Print information

        Returns: (dict) dict of alignment of the chains with the chains as keys

        '''
        self.alignment = Align.Align(self.target, self.template, output_dir=self.output_dir)

        # self.alignment = Align.Align2(target = self.target, template=self.template, output_dir=self.output_dir)
        # self.alignment.align_templates()

        if verbose:
            print('\tSuccessfully created alignment file')

    def write_ini_script(self):
        ''' Write the scipt that modeller uses for creating the initial model'''
        os.chdir(os.path.dirname(PANDORA.PANDORA_path))
        Modelling_functions.write_ini_script(self.target, self.template, self.alignment.alignment_file, self.output_dir)

    def create_initial_model(self, python_script = 'cmd_modeller_ini.py', verbose = True):
        ''' Run modeller to create the initial model. Modeller can only output files in its work directory
            (why though?), so the current work directory is changed to the output dir and later changed back the the
            old working dir.


        Args:
            python_script:  (string) path to script that performs the modeller modelling. cmd_modeller_ini.py
            verbose:  (bool) Print information

        Returns: (BIO.PDB object) self.target.initial_model

        '''
        # Change working directory
        os.chdir(self.output_dir)
        # Run Modeller
        os.popen('python %s' %python_script).read()
        # Load initial model into target object
        self.target.initial_model = PDBParser(QUIET=True).get_structure(self.target.id, self.target.id + '.ini')
        # Change working directory back
        os.chdir(os.path.dirname(PANDORA.PANDORA_path))
        if verbose:
            print('\tSuccessfully created the initital model')

    def run_modeller(self, python_script='cmd_modeller.py', benchmark=False, pickle_out=True, verbose=True, keep_IL=False):
        ''' Perform the homology modelling.

        Args:
            python_script: (string) path to script that performs the modeller modelling. cmd_modeller.py
            benchmark: (bool) Perform L-RMSD calculations? only works if the target id is an existing pdb id
            pickle_out: (bool) Save a .pkl with the results
            verbose:  (bool) Print information

        Returns: (list) list of Model objects

        '''

        if verbose:
            print('\tPerforming homology modelling of %s on %s...' %(self.target.id, self.template.id))
        t0 = time.time()
        self.results = Modelling_functions.run_modeller(self.output_dir, self.target, python_script=python_script,
                                                        benchmark=benchmark, pickle_out=pickle_out, keep_IL=keep_IL)
        if verbose:
            print('\n\tModelling was successfull and took %s seconds' %(round(time.time() - t0, 2)))

    def anchor_contacts(self, verbose=True):
        ''' Calculate anchor contacts and writes a contacts.list file that modeller uses for restraints.

        Args:
            verbose: (bool) Print information

        '''

        if verbose:
            print('\tCalculating peptide anchor residue constraints...')
        self.target.calc_anchor_contacts()
        #    Write output file
        with open(self.output_dir + '/contacts_' + self.target.id + '.list', 'w') as f:
                for i in self.target.anchor_contacts:
                    f.write('\t'.join('%s' % x for x in i) + '\n')

    def write_modeller_script(self, n_models = 20, n_jobs=None, stdev = 0.1):
        ''' Write the script that modeller uses for the final homology modelling.

        Args:
            n_models: (int) number of models that Pandora generates
            n_jobs: (int) number of parallel jobs. Is recommended to use as many jobs as the number of models: less will result in
                a slower run, more will not add any benefit but might occupy cores unnecessarily.
            stdev: (float) standard deviation of modelling restraints. Higher = more flexible restraints.

        '''
        Modelling_functions.write_modeller_script(self.target, self.template, self.alignment.alignment_file, 
                                                  self.output_dir, n_models=n_models, n_jobs=n_jobs, stdev=stdev)

    def __log(self, target_id, template_id, error, logfile = PANDORA.PANDORA_data + '/outputs/Pandora_log.txt', verbose=True):
        ''' Keeps track of what goes wrong while parsing

        Args:
            target_id: (str): ID of target structure
            template_id: (str): ID of template structure
            error: (str): error to append to log file
            logfile: (str): path to logfile
            verbose: (bool): print error?
        '''

        # Create log file
        if not os.path.exists(logfile):
            with open(logfile, 'w') as f:
                f.write('Target\tTemplate\tError\n')

        if verbose:
            print('\t' + error)
        with open(logfile, 'a') as f:
            f.write('%s\t%s\t%s\n' % (target_id, template_id, error))

    def model(self, n_models=20, n_jobs=None, stdev=0.1, seq_based_templ_selection = False, benchmark=False, verbose=True):
        ''' Wrapper function that combines all modelling steps.

        Args:
            n_models: (int) number of models modeller generates per run
            stdev: (float) standard deviation of modelling restraints. Higher = more flexible restraints.
            seq_based_templ_selection: (bool) Use template selection based on template sequences instead of allele.
            benchmark: (bool) Perform L-RMSD calculations? only works if the target id is an existing pdb id
            verbose: (bool) Print information

        Returns:

        '''

        if verbose:
            print('\nModelling %s...\n' %self.target.id)

        # Make sure we're in the root directory
        os.path.dirname(PANDORA.PANDORA_path)

        # Find the best template structure given the Target
        try:
            self.find_template(seq_based_templ_selection, benchmark=benchmark, verbose=verbose)
        except:
            self.__log(self.target.id, 'None', 'Could not find a template')
            raise Exception('Could not find a template')

        # Prepare the output directory
        try:
            self.prep_output_dir()
        except:
            self.__log(self.target.id, self.template.id, 'Failed creating output directory')
            raise Exception('Failed creating output directory')

        # Perform sequence alignment. This is used to superimpose the target on the template structure in later steps
        try:
            self.align(verbose=verbose)
        except:
            self.__log(self.target.id, self.template.id, 'Failed aligning target and template')
            raise Exception('Failed aligning target and template')

        # Prepare the scripts that run modeller
        try:
            self.write_ini_script()
        except:
            self.__log(self.target.id, self.template.id, 'Failed writing .ini script')
            raise Exception('Failed writing .ini script')

        # Run modeller to create the initial model
        try:
            self.create_initial_model(verbose=verbose)
        except:
            self.__log(self.target.id, self.template.id, 'Failed creating initial model with modeller')
            raise Exception('Failed creating initial model with modeller')

        # Calculate anchor restraints
        try:
            self.anchor_contacts(verbose=verbose)
        except:
            self.__log(self.target.id, self.template.id, 'Failed calculating anchor restraints')
            raise Exception('Failed calculating anchor restraints')

        # prepare the scripts that run modeller
        try:
            self.write_modeller_script(n_models=n_models, n_jobs=n_jobs, stdev=stdev)
        except:
            self.__log(self.target.id, self.template.id, 'Failed preparing the modeller script')
            raise Exception('Failed preparing the modeller script')

        # Do the homology modelling
        try:
            self.run_modeller(benchmark=benchmark, verbose=verbose, keep_IL=self.keep_IL)
        except:
            self.__log(self.target.id, self.template.id, 'Failed running modeller')
            raise Exception('Failed running modeller')


        if verbose and benchmark:
            try:
                print('\n\tModel\t\t\t\tMolpdf\t\tL-RMSD\t\tcore L-RMSD')
                for m in self.results:
                    try:
                        print('\t%s\t\t%s\t\t%s\t\t%s' % (
                            os.path.basename(m.model_path).replace('.pdb', ''), round(float(m.moldpf), 4),
                            round(float(m.lrmsd), 4), round(float(m.core_lrmsd), 4)))
                    except AttributeError:
                        try:
                            print('\t%s\t\t%s\t\t%s' % (
                                os.path.basename(m.model_path).replace('.pdb', ''), round(float(m.moldpf), 4),
                                round(float(m.lrmsd), 4)))
                        except AttributeError:
                            print('\t%s\t\t%s' % (
                                os.path.basename(m.model_path).replace('.pdb', ''), round(float(m.moldpf), 4)))
            except:
                self.__log(self.target.id, self.template.id, 'Could not calculate L-RMSD')
                raise Exception('Could not calculate L-RMSD')

        elif verbose and not benchmark:
            print('\n\tModel\t\t\t\tMolpdf')
            for m in self.results:
                print('\t%s\t\t%s' %(os.path.basename(m.model_path).replace('.pdb', ''), round(float(m.moldpf), 4)))

        self.__log(self.target.id, self.template.id, 'Successfully modelled %s models' %n_models)








