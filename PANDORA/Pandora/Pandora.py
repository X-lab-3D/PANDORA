
import PANDORA
from PANDORA.Pandora import Find_template
from PANDORA.Database import Database as Database
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Align
import os

class Pandora:

    def __init__(self, target, database):
        self.target = target
        self.template = None
        self.database = database
        self.output_dir = PANDORA.PANDORA_data.replace(' ', '\\ ') + '/outputs/'
        # self.alignment
        # self.ini_modeller_script
        # self.initial_model
        # self.modeller_script
        # self.results
        pass

    def find_template(self):
        self.template = Find_template.find_template(self.target, self.database)
        # create an output directory
        try:
            self.output_dir = '%s%s_%s' %(self.output_dir, self.template.PDB_id, self.target.PDB_id)
            os.system('mkdir %s' % self.output_dir)
        except:
            pass

    def find_anchors(self):
        pass

    def align(self):
        self.alignment = Align.Align(self.target, self.template)

    def write_ini_script(self):
        pass

    def run_modeller(self):
        pass

    def anchor_contacts(self):
        pass

    def write_modeller_script(self):
        pass

    def model(self):
        pass



# db = Database.Database()
# db.construct_database(MHCI=False)
#
# target = PMHC.Target('1IAK', ['MH2-AA*02', 'H2-ABk'], 'STDYGILQINSRW', MHC_class='II')
#
# mod = Pandora(target, db)
# mod.find_template()
#
# mod.target.anchors = [3,6,8,11]
# mod.template.anchors = [4,7,9,12]
#
# mod.align()
#
# mod.target.info()
# mod.alignment.alignment_file
