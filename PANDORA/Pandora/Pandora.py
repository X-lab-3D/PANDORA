
import PANDORA
from PANDORA.Pandora import Align
from PANDORA.Pandora import Modelling_functions
import time
import os
from Bio.PDB import PDBParser
from datetime import datetime


class Pandora:

    def __init__(self, target, database=None, template=None, output_dir=PANDORA.PANDORA_data + '/outputs'):
        '''__init__(self, target, database=None, template=None, output_dir=PANDORA.PANDORA_data + '/outputs')
        '''
        self.target = target
        self.template = template
        self.database = database
        self.output_dir = output_dir

        if database is None and template is None:
            raise Exception('Provide a Database object so Pandora can find the best suitable template structure for '
                            'modelling. Alternatively, you can specify a user defined Template object.')

    def find_template(self, best_n_templates=1, benchmark=False, verbose=True):
        ''' Find the best template structure given a Target object

        Args:
            best_n_templates (int, optional): how many template structures are
                used for modelling. The best n are used. Defaults to 1.
            benchmark (bool): Perform L-RMSD calculations? only works if the
                target id is an existing pdb id. Defaults to False.
            verbose: (bool): Print information. Defaults to True.

        '''

        if verbose:
            print('\tTarget MHC Class: %s' % self.target.MHC_class)
            print('\tTarget Allele:  %s' % self.target.allele_type)
            print('\tTarget Peptide: %s' % self.target.peptide)
            print('\tTarget Anchors: %s,%s\n' % (self.target.anchors[0],self.target.anchors[1]))

        if self.template is None: # Only find the best template if the user didn't specify one
            # if verbose and self.target.M_chain_seq != '' and seq_based_templ_selection:
            #     print('\tUsing sequence based template selection')
            if verbose:
                print('\tUsing allele type based template selection')
            # Find the best template. If the target already exists in the database,
            # also consider the initial loop model as a model
            self.template, self.pept_ali_scores, self.keep_IL = Modelling_functions.find_template(self.target,
                                                                    self.database,
                                                                    best_n_templates=best_n_templates,
                                                                    benchmark=benchmark)
            self.target.templates = [i.id for i in self.template]
            if verbose:
                print('\tSelected template structure (%s): %s' %(len(self.template), [i.id for i in self.template]))

        else:
            if verbose:
                if type(self.template)==list:
                    print('\tUser defined template structure (%s): %s' %(len(self.template), [i.id for i in self.template]))
                else:
                    print('\tUser defined template structure: %s' %self.template.id)
            # Check if the target structure and template structure are the same.
            if type(self.template)==list:
                self.keep_IL = any(Modelling_functions.check_target_template(self.target, tmpl) for tmpl in self.template)
            else:
                self.keep_IL = Modelling_functions.check_target_template(self.target, self.template)
            # determine peptide alignment scores of the target and the template(s)
            self.pept_ali_scores = []

            if type(self.template)==list:
                for templ in self.template:
                    if self.target.id == 'I':
                        score = Modelling_functions.score_peptide_alignment_MHCI(self.target, templ, 'PAM30')
                        self.pept_ali_scores.append((score, templ.peptide, templ.id))
                    elif self.target.id == 'II':
                        score = Modelling_functions.score_peptide_alignment_MHCII(self.target, templ, 'PAM30')
                        self.pept_ali_scores.append((score, templ.peptide, templ.id))
            else:
                if self.target.id == 'I':
                    score = Modelling_functions.score_peptide_alignment_MHCI(self.target, self.template, 'PAM30')
                    self.pept_ali_scores.append((score, self.template.peptide, self.template.id))
                elif self.target.id == 'II':
                    score = Modelling_functions.score_peptide_alignment_MHCII(self.target, self.template, 'PAM30')
                    self.pept_ali_scores.append((score, self.template.peptide, self.template.id))
            self.pept_ali_scores = self.pept_ali_scores[:best_n_templates]



        if verbose:
            if type(self.template)==list:
                print('\tTemplates Allele:  %s' %([i.allele_type for i in self.template]))
                print('\tTemplates Peptide: %s' %([i.peptide for i in self.template]))
                print('\tTemplates Anchors: %s\n' %([i.anchors for i in self.template]))
            else:
                print('\tTemplate Allele:  %s' %self.template.allele_type)
                print('\tTemplate Peptide: %s' %self.template.peptide)
                print('\tTemplate Anchors: %s\n' %self.template.anchors)

        #TODO: remove this line only after implementing issue #32.
        if type(self.template)==list:
            self.template = self.template[0]

    def prep_output_dir(self):
        ''' Create an output directory and move the template pdb there
            Uses self.output_dir (str): Path to output directory.
                Defaults to <PANDORA_location>/PANDORA_files/data/outputs.

        Args:
            None


        '''

        # create an output directory
        try:
            self.output_dir = '%s/%s_%s' %(self.output_dir, self.target.id, self.template.id)
            #self.output_dir = '%s/%s_%s' % (self.output_dir, self.target.id, '_'.join([i.id for i in self.template]))

            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
        except:
            raise Exception('A problem occurred while creating output directory')

        #for templ in self.template:

        if os.path.isfile(self.template.pdb_path):
            os.system('cp %s %s/%s.pdb' %(self.template.pdb_path, self.output_dir, self.template.id))
        else:
            print('Template object could not be found. Please check the path: %s.' %self.template.pdb_path)
            print('If the path is not available, you can use Database.repath.')
            raise Exception('Template file not found.')

        # dd/mm/YY H:M:S
        #date_time = datetime.now().strftime("%d-%m-%Y-%H-%M-%S")
        # Create the output dir of the specific case
        #self.output_dir = '%s/%s_%s_%s' %(output_dir, self.target.id, '_'.join([i.id for i in self.template]), date_time)
        #if not os.path.exists(self.output_dir):
        #    os.makedirs(self.output_dir)

        # copy the template structure to the output file
        #for t in self.template:
        os.system('cp %s %s/%s.pdb' %(self.template.pdb_path, self.output_dir, self.template.id))

    def align(self, verbose=True):
        ''' Create the alignment file for modeller.

        Args:
            verbose: (bool): Print information

        '''
        self.alignment = Align.Align(self.target, self.template, output_dir=self.output_dir)

        # self.alignment = Align.Align2(target = self.target, template=self.template, output_dir=self.output_dir)
        # self.alignment.align_templates()

        if verbose:
            print('\tSuccessfully created alignment file')

    def write_ini_script(self):
        ''' Write the python scipt that modeller uses for creating the initial model'''
        os.chdir(os.path.dirname(PANDORA.PANDORA_path))
        Modelling_functions.write_ini_script(self.target, self.template, self.alignment.alignment_file, self.output_dir)

    def create_initial_model(self, python_script = 'cmd_modeller_ini.py', verbose = True):
        ''' Run modeller to create the initial model. Modeller can only output files in its work directory
            (why though?), so the current work directory is changed to the output dir and later changed back the the
            old working dir.


        Args:
            python_script:  (str): path to script that performs the modeller modelling. Default = cmd_modeller_ini.py
            verbose:  (bool): Print information. Default = True

        '''
        # Change working directory
        os.chdir(self.output_dir)
        # Run Modeller
        os.popen('python %s > modeller_ini.log' %python_script).read()

        try:
            # Load initial model into target object
            self.target.initial_model = PDBParser(QUIET=True).get_structure(self.target.id, self.target.id + '.ini')
        except FileExistsError:
            # If the file does not exist, raise an exception to prompt the user to check MODELLER installation
            raise Exception('.ini file could not be modelled. Please check modeller_ini.log. Is your MODELLER correctly installed?')

        # Change working directory back
        os.chdir(os.path.dirname(PANDORA.PANDORA_path))
        if verbose:
            print('\tSuccessfully created the initital model')

    def run_modeller(self, python_script='cmd_modeller.py', benchmark=False, pickle_out=True, verbose=True,
                     keep_IL=False, RMSD_atoms=['C', 'CA', 'N', 'O']):
        ''' Perform the homology modelling of a target structure on template model(s). Models are saved in the output
            directory and in pandora.results[].

        Args:
            python_script: (str): path to script that performs the modeller modelling. Default = cmd_modeller.py
            benchmark: (bool): Perform L-RMSD calculations? only works if the target id is an existing pdb id.
                                Default = False
            pickle_out: (bool): Save a .pkl with the results. Default = True
            verbose: (bool): Print information. Default = True
            keep_IL: (bool): Keep the initial homology model (non optimized loops). Default = False
            RMSD_atoms: (list[str]): atoms used for the L-RMSD calculation. Default = ['C', 'CA', 'N', 'O'], which is
                                    the backbone

        '''

        if verbose:
            if type(self.template)==list:
                print('\tPerforming homology modelling of %s on %s...' %(self.target.id, '_'.join([t.id for t in self.template])))
            else:
                print('\tPerforming homology modelling of %s on %s...' %(self.target.id, self.template.id))
        t0 = time.time()
        self.results = Modelling_functions.run_modeller(self.output_dir, self.target, python_script=python_script,
                                                        benchmark=benchmark, pickle_out=pickle_out, keep_IL=keep_IL,
                                                        RMSD_atoms=RMSD_atoms)
        if verbose:
            print('\n\tModelling was successfull and took %s seconds' %(round(time.time() - t0, 2)))

    def anchor_contacts(self, verbose=True):
        ''' Calculate anchor contacts and writes a contacts.list file that modeller uses for restraints.

        Args:
            verbose: (bool): Print information. Default = True

        '''

        if verbose:
            print('\tCalculating peptide anchor residue constraints...')
        self.target.calc_anchor_contacts()
        #    Write output file
        with open(self.output_dir + '/contacts_' + self.target.id + '.list', 'w') as f:
                for i in self.target.anchor_contacts:
                    f.write('\t'.join('%s' % x for x in i) + '\n')

    def write_modeller_script(self, n_loop_models=20, n_homology_models = 1, loop_refinement='slow',
                              n_jobs=None, stdev=0.1, helix=False, sheet=False):
        ''' Write the script that modeller uses for the final homology modelling. Most modelling settings are set in
            this script.

        Args:
            n_loop_models: (int): number of loop refinement models PANDORA will generate. Default = 20
            n_homology_models: (int): number of generated homology models PANDORA generates. Default = 1
            loop_refinement: (str): levels of loop refinements. Default = slow. Supported: very_fast, fast, slow,
                                    very_slow, slow_large.
            n_jobs: (int): number of parallel jobs. Is recommended to use as many jobs as the number of models:
            less will result in a slower run, more will not add any benefit but might occupy cores unnecessarily.
            stdev: (float): standard deviation of modelling restraints that is allowed during modelling. Default = 0.1.
                A stdev of 0.2 is recommended for modelling MHC II
            helix: (bool/list): False if no alpha-helix must be modelled. Otherwise, a list of the alpha helix start
                and end-positions as integers. I.e. [3,8] for a helix between peptide residue 3 and 8.
            sheet: (bool/list): False if no beta-sheet must be modelled. Otherwise, a list containing: start position
                of B-sheet 1, start position of B-sheet 2 and the length of the B-sheet in h-bonds. For example:
                ["O:2:P","N:54:M",2] for a parallel B-sheet; The sheet starts at the Oxigen atom of the 2nd residue of
                chain P and at the Nitrogen of the 54th residue of chain M and has a length of 2 H-bonds. Or;
                ["N:6:P", "O:13:P", -3], with -3 denoting an anti-parallel B-sheet with a length of 3 H-bonds.

        '''

        Modelling_functions.write_modeller_script(self.target, self.template, self.alignment.alignment_file,
                                                  self.output_dir, n_loop_models=n_loop_models,
                                                  n_homology_models=n_homology_models, loop_refinement=loop_refinement,
                                                  n_jobs=n_jobs, stdev=stdev, helix=helix, sheet=sheet)

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

    def model(self, n_loop_models=20, n_homology_models=1,
              best_n_templates=1, n_jobs=None, loop_refinement='slow', pickle_out=False,
              stdev=0.1, benchmark=False, verbose=True, helix=False, sheet=False, RMSD_atoms=['C', 'CA', 'N', 'O']):
        '''model(self, output_dir=PANDORA.PANDORA_data + '/outputs', n_loop_models=20, n_homology_models=1, best_n_templates=1, n_jobs=None, loop_refinement='slow', pickle_out=False,stdev=0.1, benchmark=False, verbose=True, helix=False, sheet=False, RMSD_atoms=['C', 'CA', 'N', 'O'])
        Wrapper function that combines all modelling steps.

        Args:
            benchmark: (Optional, bool) If True, performs L-RMSD calculations with target strcutre.
                Only works if the target id is present in the template set.
                Defaults to False.
            helix (Optional, False or list): List of integers. Contains starting and ending
                position of a predicted alpha-helix in the peptide. Defaults to False.
            loop_refinement (Optional, str): Type of MODELLER loop refinement to apply.
                Available options are: very_fast,fast,slow,very_slow,slow_large.
                Defaults to 'slow'.
            n_loop_models (Optional, int): number of models modeller generates per run.
                Defaults to 20.
            n_homology_models (Optional, int): number of initial peptide homology models to generate.
                Defaults to 1.
            n_jobs (Optional, int or None): Number of parallel loop model jobs.
                Setting it higher than n_loop_models gives no computational time advantage.
                Recommended to change only when producing high number of loop models
                for one peptide. Defaults to None.
            output_dir (Optional, str): Path to output directory.
                Defaults to PANDORA.PANDORA_data + '/outputs'.
            pickle_out (Optional, bool): If True, saves a pickle file containing the
                PANDORA.PMHC.Model objects for the generated models in the
                output directory. Defaults to False.
            RMSD_atoms (Optional, list): list of atoms to use for final RMSD calculation.
                Works only if benchmark==True. Defaults to ['C', 'CA', 'N', 'O']
            stdev (Optional, float): standard deviation of modelling restraints.
                A higher stdev means more flexible restraints. Defaults to 0.1.
            sheet (Optional, False or list): List containing: start position of B-sheet 1,
                start position of B-sheet 2 and the length of the B-sheet in h-bonds.
                For example: ["O:2:P","N:54:M",2] for a parallel B-sheet; The sheet starts
                at the Oxigen atom of the 2nd residue of chain P and at the Nitrogen of the 54th residue of
                chain M and has a length of 2 H-bonds. Or; ["N:6:P", "O:13:P", -3], with -3 denoting an
                anti-parallel B-sheet with a length of 3 H-bonds.
            verbose (Optional, bool): If True, print modelling information. Defaults to True.

        Returns:
            None

        '''

        if verbose:
            print('\nModelling %s...\n' %self.target.id)

        # Make sure we're in the root directory
        os.path.dirname(PANDORA.PANDORA_path)

        # Find the best template structure given the Target
        try:
            self.find_template(best_n_templates=best_n_templates, benchmark=benchmark, verbose=verbose)
        except:
            self.__log(self.target.id, 'None', 'Could not find a template')
            raise Exception('Could not find a template')

        print('###############')
        print('TEMPLATE: ', self.template.id)
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
        except Exception:
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
            self.write_modeller_script(n_loop_models=n_loop_models, n_homology_models=n_homology_models,
                                       loop_refinement=loop_refinement, n_jobs=n_jobs,
                                       stdev=stdev, helix=helix, sheet=sheet)
        except:
            self.__log(self.target.id, self.template.id, 'Failed preparing the modeller script')
            raise Exception('Failed preparing the modeller script')

        # Do the homology modelling
        try:
            self.run_modeller(benchmark=benchmark, verbose=verbose, keep_IL=self.keep_IL,
                              RMSD_atoms=RMSD_atoms, pickle_out=pickle_out)
        except:
            self.__log(self.target.id, self.template.id, 'Failed running modeller')
            raise Exception('Failed running modeller')


        # if verbose and benchmark:
        #     try:
        #         print('\n\tModel\t\t\t\tMolpdf\t\tL-RMSD\t\tcore L-RMSD')
        #         molsort = sorted(self.results, key=lambda m: float(m.molpdf))
        #         for m in molsort:
        #             try:
        #                 print('\t%s\t\t%s\t\t%s\t\t%s' % (
        #                     os.path.basename(m.model_path).replace('.pdb', ''), round(float(m.molpdf), 4),
        #                     round(float(m.lrmsd), 4), round(float(m.core_lrmsd), 4)))
        #             except AttributeError:
        #                 try:
        #                     print('\t%s\t\t%s\t\t%s' % (
        #                         os.path.basename(m.model_path).replace('.pdb', ''), round(float(m.molpdf), 4),
        #                         round(float(m.lrmsd), 4)))
        #                 except AttributeError:
        #                     print('\t%s\t\t%s' % (
        #                         os.path.basename(m.model_path).replace('.pdb', ''), round(float(m.moldpf), 4)))


        #     except:
        #         self.__log(self.target.id, self.template.id, 'Could not calculate L-RMSD')
        #         raise Exception('Could not calculate L-RMSD')

        # elif verbose and not benchmark:
        if verbose:
            print('\n\tModel\t\t\t\tMolpdf')
            for m in self.results:
                print('\t%s\t\t%s' %(os.path.basename(m.model_path).replace('.pdb', ''), round(float(m.molpdf), 4)))

        if type(self.template)==list:
            self.__log(self.target.id, self.template.id, 'Successfully modelled %s models' %(n_homology_models*n_loop_models))
        else:
            self.__log(self.target.id, self.template.id, 'Successfully modelled %s models' %(n_homology_models*n_loop_models))
