
import PANDORA
from PANDORA.PMHC import PMHC
from PANDORA.Database import Database
from PANDORA.Pandora import Align
from PANDORA.Pandora import Modelling_functions
import time
import os
from Bio.PDB import PDBParser


class Pandora:

    def __init__(self, target, database = None, template = None):
        self.target = target
        self.template = template
        self.database = database
        self.output_dir = PANDORA.PANDORA_data + '/outputs/'

        if database == None and template == None:
            raise Exception('Provide a Database object so Pandora can find the best suitable template structure for '
                            'modelling. Alternatively, you can specify a user defined Template object.')

    def find_template(self, verbose = True):
        ''' Find the best template structure given a Target object '''
        if self.template == None:
            self.template = Modelling_functions.find_template(self.target, self.database)
            self.target.templates = [self.template.id]
            if verbose:
                print('\tSelected template structure: %s' %self.template.id)
        elif verbose:
                print('\tUser defined template structure: %s' %self.template.id)
        if verbose:
            print('\tTemplate Allele:  %s' % self.template.allele_type)
            print('\tTemplate Peptide: %s' % self.template.peptide)
            print('\tTemplate Anchors: %s\n' % self.template.anchors)

    def prep_output_dir(self):
        ''' Create an output directory and move the template pdb there
        '''
        # create an output directory
        try:
            self.output_dir = '%s%s_%s' %(self.output_dir, self.template.id, self.target.id)
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
        except:
            pass
        os.system('cp %s %s/%s.pdb' %(self.template.pdb_path, self.output_dir, self.template.id))

    def align(self, verbose = True):
        self.alignment = Align.Align(self.target, self.template)
        if verbose:
            print('\tSuccessfully created alignment file')

    def write_ini_script(self):
        os.chdir(os.path.dirname(PANDORA.PANDORA_path))
        Modelling_functions.write_ini_script(self.target, self.template, self.alignment.alignment_file, self.output_dir)

    def create_initial_model(self, python_script = 'cmd_modeller_ini.py', verbose = True):
        ''' Run modeller given a python script (cmd_modeller_ini.py or cmd_modeller.py). Modeller can only output files
        in its work directory (why though?), so the current work directory is changed to the output dir and later
        changed back the the old working dir.

        :param python_script: (string) path to script that performs the modeller modelling. cmd_modeller_ini.py
        :return:
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

    def run_modeller(self, python_script='cmd_modeller.py', benchmark=False, pickle_out=True, verbose = True):
        ''' Perform the homology modelling.

        :param python_script: (string) path to script that performs the modeller modelling. cmd_modeller.py
        :return:
        '''
        if verbose:
            print('\tPerforming homology modelling of %s on %s...' %(self.target.id, self.template.id))
        t0 = time.time()
        self.results = Modelling_functions.run_modeller(self.output_dir, self.target, python_script=python_script,
                                                        benchmark=benchmark, pickle_out=pickle_out)
        if verbose:
            print('\n\tModelling was successfull and took %s seconds' %(round(time.time() - t0, 2)))

    def anchor_contacts(self, verbose=True):
        """ Calculate anchor contacts"""
        if verbose:
            print('\tCalculating peptide anchor residue constraints...')
        self.target.calc_anchor_contacts()
        #    Write output file
        with open(self.output_dir + '/contacts_' + self.target.id + '.list', 'w') as f:
                for i in self.target.anchor_contacts:
                    f.write('\t'.join('%s' % x for x in i) + '\n')

    def write_modeller_script(self, n_models = 10, stdev = 0.1):
        Modelling_functions.write_modeller_script(self.target, self.template, self.alignment.alignment_file, self.output_dir, n_models=n_models, stdev=stdev)

    def model(self, n_models=10, stdev=0.1, benchmark=False, verbose=True):

        if verbose:
            print('\nModelling %s...\n' %self.target.id)
            print('\tTarget MHC Class: %s' % self.target.MHC_class)
            print('\tTarget Allele:  %s' % self.target.allele_type)
            print('\tTarget Peptide: %s' % self.target.peptide)
            print('\tTarget Anchors: %s\n' % self.target.anchors)

        # Make sure we're in the root directory
        os.path.dirname(PANDORA.PANDORA_path)
        # Find the best template structure given the Target
        self.find_template(verbose=verbose)
        # Prepare the output directory
        self.prep_output_dir()
        # Perform sequence alignment. This is used to superimpose the target on the template structure in later steps
        self.align(verbose=verbose)
        # Prepare the scripts that run modeller
        self.write_ini_script()
        # Run modeller to create the initial model
        self.create_initial_model(verbose=verbose)
        # Calculate anchor restraints
        self.anchor_contacts(verbose=verbose)
        # prepare the scripts that run modeller
        self.write_modeller_script(n_models=n_models, stdev=stdev)
        # Do the homology modelling
        self.run_modeller(benchmark=benchmark, verbose=verbose)

        if verbose and benchmark:
            print('\n\tModel\t\t\t\tMolpdf\t\tL-RMSD\t\tcore L-RMSD')
            for m in self.results:
                print('\t%s\t\t%s\t\t%s\t\t%s' % (
                os.path.basename(m.model_path).replace('.pdb', ''), round(float(m.moldpf), 4),
                round(float(m.lrmsd), 4), round(float(m.core_lrmsd), 4)))
        elif verbose and not benchmark:
            print('\n\tModel\t\t\t\tMolpdf')
            for m in self.results:
                print('\t%s\t\t%s' %(os.path.basename(m.model_path).replace('.pdb', ''), round(float(m.moldpf), 4)))

db = Database.Database()
# db.construct_database(clean=False)
# db.save('Pandora_MHCI_and_MHCII_data')
db = db.load('Pandora_MHCI_and_MHCII_data')





target = PMHC.Target('1DLH',
                     db.MHCII_data['1DLH'].allele_type,
                     db.MHCII_data['1DLH'].peptide,
                     # chain_seq = db.MHCII_data['1DLH'].chain_seq,
                     M_chain_seq = db.MHCII_data['1DLH'].M_chain_seq,
                     N_chain_seq = db.MHCII_data['1DLH'].N_chain_seq,
                     MHC_class = 'II',
                     anchors = db.MHCII_data['1DLH'].anchors)

target = PMHC.Target('1A1M',
                     db.MHCI_data['1A1M'].allele_type,
                     db.MHCI_data['1A1M'].peptide,
                     M_chain_seq = db.MHCI_data['1A1M'].M_chain_seq,
                     # MHC_class = 'II',
                     anchors = db.MHCI_data['1A1M'].anchors)


mod = Pandora(target, db)
mod.model(n_models=5, stdev=0.1, benchmark=True)





## --------------------------------------------------------------------------------




# with open('benchmark_II.csv', 'w') as f:
#     f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,L-RMSD,core L-MRSD\n')
#
# for k in db.MHCII_data:
#     try:
#         t0 = time.time()
#         print('Modelling %s' %db.MHCII_data[k].id)
#         target = PMHC.Target(db.MHCII_data[k].id, db.MHCII_data[k].allele, db.MHCII_data[k].peptide, chain_seq= db.MHCII_data[k].chain_seq, MHC_class= db.MHCII_data[k].MHC_class, anchors=db.MHCII_data[k].anchors)
#         mod = Pandora(target, db)
#         mod.model(benchmark=True)
#
#         lmrsd = round(sum([i.lrmsd for i in mod.results])/len(mod.results), 4)
#         core_lmrsd = round(sum([i.core_lrmsd for i in mod.results])/len(mod.results), 4)
#
#         print('Mean L-RMSD: %s' %(lmrsd))
#         print('Mean core L-RMSD: %s' % (core_lmrsd))
#         print('Modelling took %s seconds\n' %(time.time() - t0))
#
#         with open('benchmark_II.csv', 'a') as f:
#             f.write('%s,%s,%s,%s,%s,%s,%s,%s,\n' %(mod.target.id, mod.target.peptide, ';'.join(mod.target.allele), mod.template.id, mod.template.peptide, ';'.join(mod.template.allele), lmrsd, core_lmrsd))
#
#     except:
#         print('Something went wrong')
# #
#
#
# with open('benchmark_I.csv', 'w') as f:
#     f.write('Target_ID,Target_peptide,Target_alleles,Template_ID,Template_peptide,Template_alleles,L-RMSD,core L-MRSD\n')
#
# for k in db.MHCI_data:
#     try:
#         t0 = time.time()
#         print('Modelling %s' %db.MHCI_data[k].id)
#         target = PMHC.Target(db.MHCI_data[k].id, db.MHCI_data[k].allele, db.MHCI_data[k].peptide, chain_seq= db.MHCI_data[k].chain_seq, MHC_class= db.MHCI_data[k].MHC_class, anchors=db.MHCI_data[k].anchors)
#         mod = Pandora(target, db)
#         mod.model(benchmark=True)
#
#         lmrsd = round(sum([i.lrmsd for i in mod.results])/len(mod.results), 4)
#         core_lmrsd = round(sum([i.core_lrmsd for i in mod.results])/len(mod.results), 4)
#
#         print('Mean L-RMSD: %s' %(lmrsd))
#         print('Mean core L-RMSD: %s' % (core_lmrsd))
#         print('Modelling took %s seconds\n' %(time.time() - t0))
#
#         with open('benchmark_I.csv', 'a') as f:
#             f.write('%s,%s,%s,%s,%s,%s,%s,%s,\n' %(mod.target.id, mod.target.peptide, ';'.join(mod.target.allele), mod.template.id, mod.template.peptide, ';'.join(mod.template.allele), lmrsd, core_lmrsd))
#
#     except:
#         print('Something went wrong')
# #










