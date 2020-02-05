#!/usr/bin/env python
import modeller as M
from modeller import automodel as MA       # Load the automodel class
env = M.environ()

class MyLoop(MA.loopmodel):
    global anchor_1
    global anchor_2
    global ID
    global renum
    with open('../../data/instructions.txt', 'r') as instr_file:
        line = instr_file.readline()
        instructions = line.split(' ')
    ID = instructions[0]
    anchor_1 = int(instructions[1])
    anchor_2 = int(instructions[2])
    renum = int(instructions[3])
    #anchor_1 = 278
    #anchor_2 = 285
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids=['M', 'P'], renumber_residues=[renum,1])
        
    def select_loop_atoms(self):
        # Two residue ranges (both will be refined simultaneously)
        return M.selection(self.residue_range('%i:P' %(anchor_1 +1), '%i:P' %(anchor_2-1)))
        # note: the residue numbers and chain IDs refer to the built model, not to the template(s).

    def special_restraints(self, aln):
        rsr = self.restraints
        atoms = self.atoms
        
        # need to add more distance restraints later :)
        contact_file = open('../../data/contacts_%s.list' %ID, 'r')
        for contact_data_line in contact_file:
            contact_data = (contact_data_line.replace(' ', '')).split('\t')
            print(contact_data[3], contact_data[2], contact_data[1], contact_data[8], contact_data[7], contact_data[6], contact_data[10])
            rsr.add(M.forms.gaussian(group=M.physical.xy_distance,feature=M.features.distance(atoms['%s:%s:%s' %(contact_data[3], contact_data[2], contact_data[1])], atoms['%s:%s:%s' %(contact_data[8], contact_data[7], contact_data[6])]),mean=float(contact_data[10]), stdev=0.1))

        contact_file.close()