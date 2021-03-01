
import PANDORA
from PANDORA.Pandora import Find_template
from PANDORA.Database import Database as Database
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Align
from PANDORA.Pandora import Write_ini_script
from PANDORA.Pandora import Write_modeller_script
import os
from Bio.PDB import PDBParser
import pickle
import time

class Pandora:

    def __init__(self, target, database):
        self.target = target
        self.template = None
        self.database = database
        self.output_dir = PANDORA.PANDORA_data + '/outputs/'
        # self.alignment
        # self.ini_modeller_script
        # self.initial_model
        # self.modeller_script
        # self.results

    def find_template(self):
        ''' Find the best template structure given a Target object '''
        self.template = Find_template.find_template(self.target, self.database)


    def prep_output_dir(self):
        ''' Create an output directory and move the template pdb there
        '''
        # create an output directory
        try:
            self.output_dir = '%s%s_%s' %(self.output_dir, self.template.PDB_id, self.target.PDB_id)
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)
        except:
            pass
        os.system('cp %s %s/%s.pdb' %(self.template.pdb_path.replace(' ', '\\ '), self.output_dir.replace(' ', '\\ '), self.template.PDB_id))


    def align(self):
        self.alignment = Align.Align(self.target, self.template)

    def write_ini_script(self):
        os.chdir(os.path.dirname(PANDORA.PANDORA_path))
        Write_ini_script.write_ini_script(self.target, self.template, self.alignment.alignment_file, self.output_dir)

    def create_initial_model(self, python_script = 'cmd_modeller_ini.py'):
        ''' Run modeller given a python script (cmd_modeller_ini.py or cmd_modeller.py). Modeller can only output files
        in its work directory (why though?), so the current work directory is changed to the output dir and later
        changed back the the old working dir.

        :param python_script: (string) path to script that performs the modeller modelling. cmd_modeller_ini.py
        :return:
        '''

        # cwd = os.getcwd()
        # Change working directory
        os.chdir(self.output_dir)
        # Run Modeller
        os.popen('python %s' %python_script).read()
        # Load initial model into target object
        self.target.initial_model = PDBParser(QUIET=True).get_structure(self.target.PDB_id, self.target.PDB_id + '.ini')
        # Change working directory back
        os.chdir(os.path.dirname(PANDORA.PANDORA_path))

    def run_modeller(self, python_script = 'cmd_modeller.py'):
        ''' Perform the homology modelling.

        :param python_script: (string) path to script that performs the modeller modelling. cmd_modeller.py
        :return:
        '''
        # cwd = os.getcwd()
        # Change working directory
        os.chdir(self.output_dir)
        # run Modeller to perform homology modelling
        os.popen('python3 %s > modeller.log' %python_script).read()
        os.chdir(os.path.dirname(PANDORA.PANDORA_path))

    def anchor_contacts(self):
        """ Calculate anchor contacts"""
        self.target.calc_anchor_contacts()
        #    Write output file
        with open(self.output_dir + '/contacts_' + self.target.PDB_id + '.list', 'w') as f:
                for i in self.target.anchor_contacts:
                    f.write('\t'.join('%s' % x for x in i) + '\n')

    def write_modeller_script(self):
        Write_modeller_script.write_modeller_script(self.target, self.template, self.alignment.alignment_file, self.output_dir)

    def model(self):

        # Make sure we're in the root directory
        os.path.dirname(PANDORA.PANDORA_path)

        # Find the best template structure given the Target
        mod.find_template()
        # mod.template.anchors = [4, 7, 9, 12]
        # Prepare the output directory
        mod.prep_output_dir()
        # Perform sequence alignment. This is used to superimpose the target on the template structure in later steps
        mod.align()
        # Prepare the scripts that run modeller
        mod.write_ini_script()
        # Run modeller to create the initial model
        mod.create_initial_model()
        # Calculate anchor restraints
        mod.anchor_contacts()
        # prepare the scripts that run modeller
        mod.write_modeller_script()
        # Do the homology modelling
        mod.run_modeller()



# db = Database.Database()
# db.construct_database(MHCI=False)



# pickle.dump(db, open( "db.pkl", "wb" ) )
db = pickle.load( open( "db.pkl", "rb" ) )

#
target = PMHC.Target('1IAK', ['MH2-AA*02', 'H2-ABk'], 'STDYGILQINSRW', MHC_class='II', anchors=[3,6,8,11])
mod = Pandora(target, db)
mod.model()
#
#
# for k in db.MHCII_data:
#
#     try:
#         t0 = time.time()
#         print('Modelling %s' %db.MHCII_data[k].PDB_id)
#         target = PMHC.Target(db.MHCII_data[k].PDB_id, db.MHCII_data[k].allele, db.MHCII_data[k].peptide, MHC_class= db.MHCII_data[k].MHC_class, anchors=db.MHCII_data[k].anchors)
#         mod = Pandora(target, db)
#         mod.model()
#         print('Modelling took %s seconds\n' %(time.time() - t0))
#     except:
#         print('Something went wrong')
#












