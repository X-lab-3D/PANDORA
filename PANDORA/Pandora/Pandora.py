
import PANDORA
from PANDORA.Pandora import Find_template
from PANDORA.Database import Database as Database
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Align
from PANDORA.Pandora import Write_ini_script
import os
from Bio.PDB import PDBParser

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
        pass

    def find_template(self):
        self.template = Find_template.find_template(self.target, self.database)


    def prep_output_dir(self):
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
        Write_ini_script.write_ini_script(self.target, self.template, self.alignment.alignment_file, self.output_dir)

    def create_initial_model(self, python_script = 'cmd_modeller_ini.py'):
        ''' Run modeller given a python script (cmd_modeller_ini.py or cmd_modeller.py). Modeller can only output files
        in its work directory (why though?), so the current work directory is changed to the output dir and later
        changed back the the old working dir.

        :param python_script:
        :return:
        '''

        cwd = os.getcwd()
        # Change working directory
        os.chdir(self.output_dir)
        # Run Modeller
        os.popen('python %s' %python_script).read()
        # Load initial model into target object
        self.target.initial_model = PDBParser(QUIET=True).get_structure(self.target.PDB_id, self.target.PDB_id + '.ini')
        # Change working directory back
        os.chdir(cwd)

    def run_modeller(self):
        pass

    def anchor_contacts(self):
        """ Calculate anchor contacts"""
        self.target.calc_anchor_contacts()
        #    Write output file
        with open(self.output_dir + '/' + self.target.PDB_id + '.list', 'w') as f:
                for i in self.target.anchor_contacts:
                    f.write('\t'.join('%s' % x for x in i) + '\n')

    def write_modeller_script(self):
        pass

    def model(self):
        pass



db = Database.Database()
db.construct_database(MHCI=False)

target = PMHC.Target('1IAK', ['MH2-AA*02', 'H2-ABk'], 'STDYGILQINSRW', MHC_class='II')

mod = Pandora(target, db)

mod.find_template()
mod.target.anchors = [3,6,8,11]
mod.template.anchors = [4,7,9,12]

mod.prep_output_dir()
mod.align()

mod.write_ini_script()
mod.create_initial_model()

mod.anchor_contacts()


for i in mod.target.anchor_contacts:
    print(i)

mod.output_dir











