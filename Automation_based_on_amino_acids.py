import os
from itertools import permutations

def obtain_seperate_rtf(ptm):
    info = ""
    apply_patch = True
    residue_topology_files = [os.path.join(os.getcwd(), "Toppar/toppar_all36_prot_modify_res.str"), 
        os.path.join(os.getcwd(), "Toppar/toppar_all36_prot_na_combined.str")]
    # with open("/home/shahielm/toppar/top_all36_prot.rtf", "r") as file:
    for topology_file in residue_topology_files:
        print(topology_file)
        with open(topology_file, "r") as file:
            data = file.readlines()
        file.close()
        res = [i for i in data if f"RESI {ptm}" in i]
        if res != "":
            part_of_ptm = False
            for line in data:
                if line.strip():
                    if "DOUBLE" in line:
                        line = line.replace("DOUBLE", "BOND  ")
                    if "DOUB" in line:
                        line = line.replace("DOUB", "BOND")
                    # Removes unnecessary information that is mentioned after an "!"
                    if "!" in line:
                        line = line.split("!", 1)[0]
                        # Some lines only have text after the "!", this strip is done to prevent a bunch of empty lines in the corrected PTM topologies
                        line = line.strip()
                        # This is true if the line has text before the "!"
                        if line != "":
                            # add the newline back which was removed with the split command. This prevents that two lines will be put on the same line
                            line += "\n"
                    if f"RESI {ptm}" in line:
                        part_of_ptm = True
                        info = line
                    elif "RESI" in line or "PRES" in line:
                        part_of_ptm = False
                    elif part_of_ptm == True and "GROUP" in line:
                        pass
                    elif part_of_ptm == True and "CMAP" in line:
                        pass
                    elif part_of_ptm == True and "DONOR" in line:
                        pass
                    elif part_of_ptm == True and "ACCEPTOR" in line:
                        pass
                    elif part_of_ptm == True:
                        info += line       
            with open(f"{ptm}.txt", "w") as output:
                output.write(info + "PATC FIRS NTER LAST CTER\n")
            output.close()

def create_edited_file(ptm):
    dictionary = {
    "HB1" : "HB ", # not in cgenff, Backbone hydrogen
    "CG321" : "CT2  ", # alphatic proton, CH2
    "HGA2" : "HA  ", # alphatic proton, CH2
    "NG2S1" : "NH1  ", # peptide nitrogen (CO=NHR)
    "HGP1" : "H   ", # polar H
    "CG2O6" : "CC   ", # carbonyl C: urea, carbonate --> urea and carbonate not mentioned solv.lib and radii.lib
    "OG2D1" : "O    ", # carbonyl O: amides, esters, [neutral] carboxylic acids, aldehydes, urea
    "NG2S2" : "NH2  ", # terminal amide nitrogen (CO=NH2)
    "HA2" : "HA ", # not in cgenff, CH2
    "HA3" : "HA ", # not in cgenff, CH3
    "CG324" : "CT2  ", # aliphatic C for CH2, adjacent to positive N (piperidine) (+) --> Could not find an exact match in the sol.lib
    "NG2P1" : "NC2  ", # N for protonated imine/Schiff's base (C=N(+)H-R, acyclic amidinium, guanidinium) --> Could not find an exact match in the solv.lib
    "HGP2" : "H   ", # polar H, +ve charge
    "CG2N1" : "C    ", # conjugated C in guanidine/guanidinium
    "CG334" : "CT3  ", # aliphatic C for methyl group (-CH3), adjacent to positive N (PROT NTER) (+) --> Could not find an exact match in the solv.lib
    "HGA3" : "HA  ", # alphatic proton, CH3
    "NG3P2" : "NH1  ", # secondary NH2+ (proline)
    "SG301" : "S    ", # sulfur C-S-S-C type
    "SG311" : "S    ", # sulphur, SH, -S-
    "HGP3" : "H   " # polar H, thiol
}

    # "CG2R61" : "CA    ", # 6-mem aromatic C
    # "HGR62" : "HA   ", # nonpolar H, neutral 6-mem planar ring C adjacent to heteroatom
    # "HGR61" : "HP   ", # aromatic H
    # "CG2R66" : "CA    " #6-mem aromatic carbon bound to F
    # #"FGR1" : "" # aromatic fluorine

# Note 1 Acetylated lysine in the toppar file is requires very little editing. Why some of the amino acids in the toppar file don't have any cgennff atom type names, I don't know.
# Issue 1: The bond in the description of NG2S1 is not present in the residue (CO=NH2). In the cgenff file. this is correct.
# Issue 2: NG2S2 is not at near one of the terminals, yet the description says it is, same explanations as issue 1.
# Issue 3: the description for CG324 does not match the structure of NMM --> The N next to this C is not the N that is protonated.
# Issue 4: the description for NG2P1 does not match the structure of NMM --> Only one of the 3 N's in NMM should have the atom type CG2N1 (because only one is protonated).
# Issue 5: the description for CG334 does not match the structure of NMM --> The atom next to this C is not the protonated N.
# Issue 6: In cgenff, NG3P2 is defined as "secondary NH2+ (proline)". In toppar file this is a protonated nitrogen. This means that the toppar file with the modified amino acids is wrong regarding this atom type. So technically this dictionary key-value pair is wrong because the modified amino acid file is wrong, so if you use a correct rtf file this might be incorrectly defined.
    
    with open(f"{ptm}.txt", "r") as file:
        data = file.readlines()
    file.close
    index = -1
    for line in data:
        index += 1
        if "ATOM" in line:
            difference = 0
            split_line = line.split()
            if split_line[2] in dictionary:
                if len(split_line[2]) < len(dictionary[split_line[2]]):
                    difference = (len(dictionary[split_line[2]]) - len(split_line[2]))
                #print(split_line[2])
                #print(dictionary[split_line[2]])
                if difference > 0:
                    data[index] = (line.replace(split_line[2] + difference*" ", dictionary[split_line[2]]))
                else:
                    data[index] = (line.replace(split_line[2], dictionary[split_line[2]]))

            # P2 (pyrophosphate phosphorus) and ON2B (Nucleic acid phosphate ester oxygen, tyr-phosphate) might also need to be included, since they are part of a different type of phosphate as seen in SEP
            elif split_line[2] in ["ON2", "P", "ON3", "ON4", "HN4"]:
                print(f"{split_line[2]} is a nucleic acid atom type which is part of a phosphate. Which is odd since SEP not being a nucleic acid.")
            # Used to easily find the atom types in the PTM topologies that are not implemented in the dictionary.  
            elif split_line[2] not in ["NH1", "H", "CT1", "HB1", "C", "O", "CT2", "CT3"]:
               print(f"{split_line[2]} not implemented in dictionary.")

    # for sentence in data:
    #     print(sentence)
    with open(f"{ptm}_edited.txt", "w") as edited_file:
        for line2 in data:
            edited_file.write(line2)
    edited_file.close()

def write_to_top_heav_lib(ptm, last_ptm, pathway):
    # Opening our text file in write only mode to write the replaced content
    with open(os.path.join(os.path.dirname(os.getcwd()), 'miniconda3/lib/modeller-10.2/modlib/top_heav.lib'), 'r') as top_lib:
        top_lib_data = top_lib.read()
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
            # print(ptm)
            # print(last_ptm)
            if ptm == last_ptm:
                # print(data[:-1])
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
            # (custom made top.lib). The second newline is to create a gap between the different topologies
            writefile.write("\n\n")
            writefile.write(data)
            print("Text replaced")
            writefile.close()

def write_to_restyp_lib(ptms, pathway):
    restyp_read_file = open(os.path.join(os.path.dirname(os.getcwd()), '/home/shahielm/miniconda3/lib/modeller-10.2/modlib/restyp.lib'), 'r')
    #restyp_read_file = open('/home/shahielm/miniconda3/lib/modeller-10.2/modlib/restyp.lib', 'r')
    restyp_lines = restyp_read_file.read()
    restyp_read_file.close()
    similar_dict = {"CIR" : ["R","Citrulline"],
                "ALY" : ["K","Acetylated Lysine"],
                "NMM" : ["R","Methylated Arginine"],
                "MLZ" : ["K","Methylated Lysine"],
                "SEP" : ["S","Phosphoserine"]
}
    counter = 0
    symbols = ["8","9", "+", "=", "?", "`", "*", "q", "U", "O", "J", "\\", ";", "'", 
            ",", "~", "!", "%", "^", "(", ")", "_", "{", "}", "|", ":", '"',"<",
            ">"]
    restyp_ptms = ""

    # Make the restyp.lib lines for the PTMs
    for ptm in ptms:
        description = similar_dict[ptm][1]
        if ptm in similar_dict.keys():
            similar = similar_dict[ptm][0]
            line = f"HETATM | {ptm}                 | {symbols[counter]} | {similar} | {ptm}  | {description}"
        else:
            line = f"HETATM | {ptm}                 | {symbols[counter]} |  | {ptm}  | {description}"
        restyp_ptms += line + "\n"
        counter += 1

    #removes empty last line 
    restyp_ptms = restyp_ptms[:-1]
    with open(pathway, 'w') as writefile:
        # If it is the last ptm that is appended to the file, then the last 
        # character will be removed. This is done since the last character 
        # is a newline character which causes an unneeded empty line.
        writefile.write(restyp_lines)
        writefile.write("\n")
        writefile.write(restyp_ptms)
        writefile.close()

ptms = ["CIR", "ALY", "NMM", "MLZ", "SEP"]
#ptms = ["F2F"]

pathway_top_heav_lib = os.path.join(os.getcwd(), "combined_top_heav_lib_file.lib")
pathway_restyp_lib = os.path.join(os.getcwd(), "combined_restyp_file.lib")

# pathway_top_heav_lib = '/home/shahielm/PANDORA/combined_top_heav_lib_file.lib'
# pathway_restyp_lib = '/home/shahielm/PANDORA/combined_restyp_file.lib'

if os.path.exists(pathway_top_heav_lib):
    os.remove(pathway_top_heav_lib)
for ptm in ptms:
    print(ptm)
    obtain_seperate_rtf(ptm)
    create_edited_file(ptm)
    write_to_top_heav_lib(ptm, ptms[-1], pathway_top_heav_lib)
write_to_restyp_lib(ptms, pathway_restyp_lib)