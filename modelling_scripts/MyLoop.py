#!/usr/bin/env python
import modeller as M
from modeller import automodel as MA       # Load the automodel class
env = M.environ()

class MyLoop(MA.loopmodel):
    global anchor_1
    global anchor_2
    anchors = []
    with open('data/instructions.txt', 'r') as instr_file:
        for i, line in enumerate(instr_file):
            anchors.append(int(line.rstrip()))
    anchor_1 = anchors[0]
    anchor_2 = anchors[1]
    #anchor_1 = 278
    #anchor_2 = 285
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['M', 'P'])#, renumber_residues=[1,1])
        self.write(file = 'rename_try.pdb')
        
    def select_loop_atoms(self):
        # Two residue ranges (both will be refined simultaneously)
        return M.selection(self.residue_range('%i:P' %anchor_1, '%i:P' %anchor_2))
        # note: the residue numbers and chain IDs refer to the built model, not to the template(s).

    def special_restraints(self, aln):
        rsr = self.restraints
        atoms = self.atoms
        
        # need to add more distance restraints later :)
        contact_file = open('data/contacts_P%i_P%i.list' %(anchor_1, anchor_2), 'r')
        for contact_data_line in contact_file:
            contact_data = (contact_data_line.replace(' ', '')).split('\t')
            print(contact_data[3], contact_data[2], contact_data[1], contact_data[8], contact_data[7], contact_data[6], contact_data[10])
            rsr.add(M.forms.gaussian(group=M.physical.xy_distance,feature=M.features.distance(atoms['%s:%s:%s' %(contact_data[3], contact_data[2], contact_data[1])], atoms['%s:%s:%s' %(contact_data[8], contact_data[7], contact_data[6])]),mean=float(contact_data[10]), stdev=0.1))

        contact_file.close()
