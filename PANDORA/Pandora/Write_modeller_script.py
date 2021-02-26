import os
import PANDORA

def write_modeller_script(target, template, alignment_file, output_dir):
    ''' Write script that refines the loops of the peptide

    :param target:
    :param template:
    :param alignment_file:
    :param output_dir:
    :return:  writes two files; Myloop.py and cmd_modeller.py that contain the info to run modeller
    '''


    anch = target.anchors

    with open(output_dir.replace('\\ ', ' ') + '/MyLoop.py', 'w') as myloopscript:
        MyL_temp = open(PANDORA.PANDORA_path + '/Pandora/MyLoop_template.py', 'r')
        for line in MyL_temp:
            if 'self.residue_range' in line:
                myloopscript.write(line % (1, anch[0]))  # write the first anchor
                for i in range(len(anch) - 1):  # Write all the inbetween acnhors if they are there
                    myloopscript.write(line % (anch[i] + 2, anch[i + 1]))
                myloopscript.write(line % (anch[-1] + 2, len(target.peptide)))  # Write the last anchor
            elif 'contact_file = open' in line:
                myloopscript.write(line % target.PDB_id)
            else:
                myloopscript.write(line)
        MyL_temp.close()


    with open(output_dir.replace('\\ ', ' ') + '/cmd_modeller.py', 'w') as modscript:
        cmd_m_temp = open(PANDORA.PANDORA_path + '/Pandora/cmd_modeller_template.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line %(os.path.basename(alignment_file)))
            elif 'knowns' in line:
                modscript.write(line %(template.PDB_id, target.PDB_id))
            else:
                modscript.write(line)
        cmd_m_temp.close()



