def write_ini_script(target, template, alignment_file, output_dir):
    ''' Writes the MyLoop.py and cmd_modeller_ini.py files. This function takes two template python scripts and fills
    in the required information: Anchor positions for the MyLoop file and structure name + alignment file for the
    cmd_modeller_ini file.

    :param target:
    :param template:
    :param alignment_file:
    :param output_dir:
    :return:
    '''

    anch = target.anchors

    with open(output_dir.replace('\\ ',' ') + '/MyLoop.py', 'w') as myloopscript:
        MyL_temp = open('PANDORA/Pandora/MyLoop_template.py', 'r')
        for line in MyL_temp:
            if 'self.residue_range' in line:
                myloopscript.write(line % (1, anch[0])) # write the first anchor
                for i in range(len(anch)-1): # Write all the inbetween acnhors if they are there
                    myloopscript.write(line % (anch[i] + 2, anch[i+1]))
                myloopscript.write(line % (anch[-1] + 2, len(template.peptide))) # Write the last anchor
            elif 'SPECIAL_RESTRAINTS_BREAK' in line:
                break
            elif 'contact_file = open' in line:
                myloopscript.write(line % template_ID)
            else:
                myloopscript.write(line)
        MyL_temp.close()

    with open(output_dir.replace('\\ ', ' ') + '/cmd_modeller_ini.py', 'w') as modscript:
        cmd_m_temp = open('PANDORA/Pandora/cmd_modeller_ini.py', 'r')
        for line in cmd_m_temp:
            if 'alnfile' in line:
                modscript.write(line % alignment_file)
            elif 'knowns' in line:
                modscript.write(line % (template.PDB_id, target.PDB_id))
            else:
                modscript.write(line)
        cmd_m_temp.close()

