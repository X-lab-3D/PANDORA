
from modeller import *
from modeller.automodel import *             # Load the automodel class
env = environ()

class MyLoop(loopmodel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['A', 'C'],
                             renumber_residues=[1,1])
    def select_loop_atoms(self):
        # Two residue ranges (both will be refined simultaneously)
        return selection(self.residue_range('3:C', '8:C'))
        # note: the residue numbers and chain IDs refer to the built model, not to the template(s).

    def special_restraints(self, aln):
        rsr = self.restraints
        atoms = self.atoms
        
        # need to add more distance restraints later :)
        contact_file = open('contacts_P2_P9.list', 'r')
        for contact_data_line in contact_file:
            contact_data = contact_data_line.split('\t')
            rsr.add(forms.gaussian(group=physical.xy_distance,feature=features.distance(atoms['%s:%s:%s' %(contact_data[3], contact_data[2], contact_data[1])], atoms['%s:%s:%s' %(contact_data[8], contact_data[7], contact_data[6])]),mean=float(contact_data[10]), stdev=0.1))

        contact_file.close()
