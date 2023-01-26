import os

def check_ptm(ptm, residue_topology_files):
    '''Checks if the given PTMs have a topology in the residue topology files.

    Args:
        ptm (str): 3 letter code for the post-translational modification.
        residue_topology_files (list): list that contains the paths to the files
            in which the PTM topologies are stored.

    Returns:
        None

        '''
    found_ptm_topology = False
    for topology_file in residue_topology_files:
        with open((topology_file), 'r') as check_top_heav_lib:
            if f"RESI {ptm}" in check_top_heav_lib.read():
                found_ptm_topology = True
    # Raise excepettion if the PTM topology was not found in any of the residue
    # topology files.
    if found_ptm_topology == False:
        raise Exception(f"{ptm} topology was not found.")

def obtain_seperate_rtf(ptm, residue_topology_files):
    '''Retrieves the PTM topologies, modifies them to get rid of unnecessary 
    information and writes these modified PTM topologies into files.

    Args:
        ptm (str): 3 letter code for the post-translational modification.
        residue_topology_files (list): list that contains the paths to the files
            in which the PTM topologies are stored.

    Returns:
        None

    '''
    
    info = ""
    # Retrieves the topologies for each of the provided PTMs from the residue
    # topology files.
    for topology_file in residue_topology_files:
        with open(topology_file, "r") as file:
            data = file.readlines()
        file.close()
        res = [i for i in data if f"RESI {ptm}" in i]
        if res != "":
            part_of_ptm = False
            for line in data:
                if line.strip():
                    # Modeller uses BOND as a keyword for any type of covalent
                    # bond. Meaning that the DOUBLE keyword needs to be replaced
                    if "DOUBLE" in line:
                        line = line.replace("DOUBLE", "BOND  ")
                    # Some PTMs have DOUBLE written as DOUB, so an extra if 
                    # statement is needed.
                    if "DOUB" in line:
                        line = line.replace("DOUB", "BOND")
                    # TRIPLE keyword needs to be replaced for the same reason as
                    # DOUBLE.
                    if "TRIPLE" in line:
                        line = line.replace("TRIPLE", "BOND")
                    # Removes unnecessary information that is mentioned after an 
                    # "!"
                    if "!" in line:
                        line = line.split("!", 1)[0]
                        # Some lines only have text after the "!", this strip 
                        # is done to prevent a bunch of empty lines in the 
                        # corrected PTM topologies
                        line = line.strip()
                        # This is true if the line has text before the "!"
                        if line != "":
                            # add the newline back which was removed with the 
                            # split command. This prevents that two lines will 
                            # be put on the same line
                            line += "\n"
                    
                    unnecessary_keywords = ["GROUP", "CMAP", "DONOR", "ACCEPTOR"]
                    # Starts storing the lines when the needed PTM topology 
                    # occurs. 
                    if f"RESI {ptm}" in line:
                        part_of_ptm = True
                        info = line
                    # Stops storing the lines when the next PTM topology occurs.
                    elif "RESI" in line or "PRES" in line:
                        part_of_ptm = False
                    # Don't store the lines that have unnecessary information
                    # (GROUP, CMAP, DONOR, ACCEPTOR).
                    elif part_of_ptm == True and any(keyword in line for keyword in unnecessary_keywords):
                        pass
                    # Keeps storing the information as long as the needed PTM 
                    # topology hasn't ended yet, while skipping the unnecessary
                    # information.
                    elif part_of_ptm == True:
                        info += line       
            # Writes the stored PTM topology along with the patches that need to be
            # applied when the residue occurs in the first or last position of
            # a chain into a file.
            with open(f"{ptm}.txt", "w") as output:
                output.write(info + "PATC FIRS NTER LAST CTER\n")
            output.close()

def create_edited_file(ptm):
    '''Checks if the given PTMs have a topology in the residue topology files.

    Args:
        ptm (str): 3 letter code for the post-translational modification.
        residue_topology_files (list): list that contains the paths to the files
            in which the PTM topologies are stored.

    Returns:
        None

        '''

    # Dictionary which contains the CHARMM36 atom type names (descriptions from 
    # /PANDORA/PTM_implementation/Toppar/top_all36_cgenff.rtf or 
    # /PANDORA/PTM_implementation/Toppar/toppar_all36_prot_na_combined.str) as 
    # the key and the CHARMM22 atom type names (desc from 
    # /miniconda3/lib/modeller-10.2/modlib/radii.lib or 
    # /miniconda3/lib/modeller-10.2/modlib/solv.lib) as the value. The spaces 
    # are needed to make sure all atom type names line up in the PTM topology. 
    dictionary = {
    "HB1" : "HB ",     # not in cgenff
    "CG321" : "CT2  ", # aliphatic C for CH2
    "HGA2" : "HA  ",   # alphatic proton, CH2
    "NG2S1" : "NH1  ", # peptide nitrogen (CO=NHR)
    "HGP1" : "H   ",   # polar H
    "CG2O6" : "CC   ", # carbonyl C: urea, carbonate --> urea and carbonate not mentioned solv.lib and radii.lib
    "OG2D1" : "O    ", # carbonyl O: amides, esters, [neutral] carboxylic acids, aldehydes, urea
    "NG2S2" : "NH2  ", # terminal amide nitrogen (CO=NH2)
    "HA2" : "HA ",     # not in cgenff
    "HA3" : "HA ",     # not in cgenff
    "CG324" : "CT2  ", # aliphatic C for CH2, adjacent to positive N (piperidine) (+) --> Could not find an exact match in the sol.lib
    "NG2P1" : "NC2  ", # N for protonated imine/Schiff's base (C=N(+)H-R, acyclic amidinium, guanidinium) --> Could not find an exact match in the solv.lib
    "HGP2" : "H   ",   # polar H, +ve charge
    "CG2N1" : "C    ", # conjugated C in guanidine/guanidinium
    "CG334" : "CT3  ", # aliphatic C for methyl group (-CH3), adjacent to positive N (PROT NTER) (+) --> Could not find an exact match in the solv.lib
    "HGA3" : "HA  ",   # alphatic proton, CH3
    "NG3P2" : "NH1  "  # secondary NH2+ (proline)
}

    with open(f"{ptm}.txt", "r") as file:
        data = file.readlines()
    file.close
    index = -1
    for line in data:
        index += 1
        if "ATOM" in line:
            difference = 0
            split_line = line.split()
            # if atom type is already implemented, correct its name formatting
            if split_line[2] in dictionary:
                if len(split_line[2]) < len(dictionary[split_line[2]]):
                    difference = (len(dictionary[split_line[2]]) - len(split_line[2]))
                if difference > 0:
                    data[index] = (line.replace(split_line[2] + difference*" ", dictionary[split_line[2]]))
                else:
                    data[index] = (line.replace(split_line[2], dictionary[split_line[2]]))
###
            # TODO: P2 (pyrophosphate phosphorus) and ON2B (Nucleic acid phosphate ester oxygen, tyr-phosphate) might also need to be included, since they are part of a different type of phosphate as seen in SEP
            # SEP
            #if the atom belongs to a nucleic acid in CHARMM22
            elif split_line[2] in ["ON2", "P", "ON3", "ON4", "HN4"]:
                print(f"WARNING: {split_line[2]} is a nucleic acid atom type which is part of a phosphate. SEP is not a nucleic acid. Please double-check your input.")
            # Used to easily find the atom types in the PTM topologies that are 
            # not implemented in the dictionary.  
            elif split_line[2] not in ["NH1", "H", "CT1", "HB1", "C", "O", "CT2", "CT3"]:
                print(f"WARNING: {split_line[2]} is a new atom type which is not implemented in the dictionary.")

    # Write down edits in a new file
    with open(f"{ptm}_edited.txt", "w") as edited_file:
        for line2 in data:
            edited_file.write(line2)
    edited_file.close()

def write_to_top_heav_lib(ptm, last_ptm, pathway):
    '''Checks if the given PTMs have a topology in the residue topology files.

    Args:
        ptm (str): 3 letter code for the post-translational modification.
        residue_topology_files (list): list that contains the paths to the files
            in which the PTM topologies are stored.

    Returns:
        None

    '''

    # Opening the text file in write only mode to write the replaced content
    with open(os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), 'miniconda3/lib/modeller-10.2/modlib/top_heav.lib'), 'r') as top_lib:
        top_lib_data = top_lib.read()
    top_lib.close()
    with open(f"{ptm}_edited.txt", "r") as ptm_file:
        data = ptm_file.read()
    ptm_file.close()
    if os.path.exists(pathway):
        with open(pathway, 'a') as appendfile:
            # Writing the replaced data in our text file
            # The newline is to create a gap between the different topologies
            appendfile.write("\n")
            # If it is the last ptm that is appended to the file, then the last 
            # character will be removed. This is done since the last character 
            # is a newline character which causes an unneeded empty line.
            if ptm == last_ptm:
                data = data[:-1]
            appendfile.write(data)
            print("Text replaced")
            appendfile.close()
    else:
        with open(pathway, 'w') as writefile:
            # Start with writing the base top_heav.lib topologies to this file
            writefile.write(top_lib_data)
            # Writing the replaced data into our text file
            # two newline characters: the first is needed so that the first line 
            # of this file will not be written on the same line as the last line 
            # in the file which this file will be appended to 
            # (custom made top.lib). The second newline is to create a gap 
            # between the different topologies.
            writefile.write("\n\n")
            writefile.write(data)
            print("Text replaced")
            writefile.close()

def write_to_restyp_lib(ptms, pathway):
    '''Checks if the given PTMs have a topology in the residue topology files.

    Args:
        ptm (str): 3 letter code for the post-translational modification.
        residue_topology_files (list): list that contains the paths to the files
            in which the PTM topologies are stored.

    Returns:
        None

    '''

    restyp_read_file = open(os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), "miniconda3/lib/modeller-10.2/modlib/restyp.lib"), 'r')
    restyp_lines = restyp_read_file.read()
    restyp_read_file.close()
    # Three-letter : [Similar amino acid, Description]
    similar_dict = {"CIR" : ["R","Citrulline"],
                "ALY" : ["K","Acetylated Lysine"],
                "NMM" : ["R","Methylated Arginine"],
                "MLZ" : ["K","Methylated Lysine"],
                "SEP" : ["S","Phosphoserine"]
                }
    counter = 0
    # Symbols that are available, but not used: *, \\, ', ", q, U, O, J
    symbols = ["8","9", "+", "=", "?", "`", ";", 
            ",", "~", "!", "%", "^", "(", ")", "_", "{", "}", "|", ":","<",
            ">"]
    restyp_ptms = ""

    # Make the restyp.lib lines for the PTMs
    for ptm in ptms:
        description = similar_dict[ptm][1]
        ##
        if ptm in similar_dict.keys():
            similar = similar_dict[ptm][0]
            line = f"HETATM | {ptm}                 | {symbols[counter]} | {similar} | {ptm}  | {description}"
        else:
            line = f"HETATM | {ptm}                 | {symbols[counter]} |  | {ptm}  | {description}"
        restyp_ptms += line + "\n"
        counter += 1

    # Removes empty last line 
    restyp_ptms = restyp_ptms[:-1]
    with open(pathway, 'w') as writefile:
        # If it is the last ptm that is appended to the file, then the last 
        # character will be removed. This is done since the last character 
        # is a newline character which causes an unneeded empty line.
        writefile.write(restyp_lines)
        writefile.write("\n")
        writefile.write(restyp_ptms)
        writefile.close()

def main():
    # Provides the PTMs (as in toppar_all36_prot_na_combined.str and toppar_all36_prot_na_combined.str) 
    # and the paths to the files needed for the retrieval of the topologies.
    ptms = ["CIR", "ALY", "NMM", "MLZ", "SEP"]
    pathway_custom_top_heav_lib = os.path.join(os.getcwd(), "combined_top_heav_lib_file.lib")
    pathway_restyp_lib = os.path.join(os.getcwd(), "combined_restyp_file.lib")
    residue_topology_files = [os.path.join(os.getcwd(), "Toppar/toppar_all36_prot_modify_res.str"), os.path.join(os.getcwd(), "Toppar/toppar_all36_prot_na_combined.str")]

    # Removes the custom top_heav_lib file so a new one can be made.
    if os.path.exists(pathway_custom_top_heav_lib):
        os.remove(pathway_custom_top_heav_lib)

    # Loop that will retrieve, modify and implement the PTMs into Modeller.
    for ptm in ptms:
        check_ptm(ptm, residue_topology_files)
        obtain_seperate_rtf(ptm, residue_topology_files)
        create_edited_file(ptm)
        write_to_top_heav_lib(ptm, ptms[-1], pathway_custom_top_heav_lib)
    write_to_restyp_lib(ptms, pathway_restyp_lib)

if __name__=='__main__':
    main()