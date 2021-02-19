
import os
import gzip
import shutil

import PANDORA
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO


def get_chainid_alleles_MHCII(pdbf):
    '''
    Takes as input an IMGT preprocessed PDB file of p:MHC II.
    Returns a dictionary containing alleles andrelative identity scores for each
    G-domain in the given pdb from the REMARK.

    Args:
        pdbf(str) : path to IMGT pdb file
    '''
    # test: multiple chains 3GJG, multiple alleles 1AO7
    ### Parsing file and extracting remarks
    with open(pdbf) as infile:
        remarks = []
        for line in infile:
            if line.startswith('REMARK 410'):
                row = [x for x in line.rstrip().split(' ') if x != '']
                del row[:2]
                remarks.append(row)
    remarks = [x for x in remarks if x != []]

    ### Dividing each remark section into a chains dictionary
    chains = {}
    flag = False
    for row in remarks:
        if row[0] == 'Chain' and row[1] == 'ID' and len(row) == 4:
            chainID = row[2][-1]
            chains[chainID] = []
            chains[chainID].append(row)
            flag = True
        elif flag == True:
            chains[chainID].append(row)

    ### Extracting MHC II Alpha and Beta chains
    mhc_a = {}  # MHC II Alpha
    mhc_b = {}  # MHC II Beta
    for chain in chains:
        try:
            if chains[chain][1][3] == 'II-ALPHA':
                mhc_a[chain] = chains[chain]
            elif chains[chain][1][3] == 'II-BETA':  # or chains[chain][1][3] == 'II-BETA':
                mhc_b[chain] = chains[chain]
        except:
            pass

    ### Extracting alleles info
    mhc_a_alleles = {x: [] for x in mhc_a}
    mhc_b_alleles = {x: [] for x in mhc_b}

    for chain in mhc_a:
        key = False
        for row in mhc_a[chain]:
            if row[0] == 'G-DOMAIN':
                try:
                    if row[3] == 'description' and row[4] == 'G-ALPHA':
                        key = 'G-ALPHA'
                    elif key:
                        if row[2] == 'gene' and row[3] == 'and' and row[4] == 'allele':
                            mhc_a_alleles[chain] += row[5:]
                        else:
                            key = False
                except IndexError:
                    pass

    for chain in mhc_b:
        key = False
        for row in mhc_b[chain]:
            if row[0] == 'G-DOMAIN':
                try:
                    if row[3] == 'description' and row[4] == 'G-BETA':
                        key = 'G-BETA'
                    elif key:
                        if row[2] == 'gene' and row[3] == 'and' and row[4] == 'allele':
                            mhc_b_alleles[chain] += row[5:]
                        else:
                            key = False
                except IndexError:
                    pass

    '''
    mhc_a_alleles_percs = {}
    for chain in mhc_a_alleles:
        mhc_a_alleles_percs[chain] = {}
        for key in mhc_a_alleles[chain]:
            mhc_a_alleles_percs[chain][key] = {}
            ### Allele info are always given with four elements: Gender, Spieces, Allele, Percentage
            for block in range(int(len(mhc_a_alleles[chain][key])/4)):
                allele = mhc_a_alleles[chain][key][2+(4*block)]
                perc = float(mhc_a_alleles[chain][key][3+(4*block)].replace('(', '').replace('%)', '').replace(',',''))
                mhc_a_alleles_percs[chain][key][allele] = perc
    '''
    return {'Alpha': mhc_a_alleles, 'Beta': mhc_b_alleles}

def unzip_pdb(ID, indir, outdir):
    """ Unzips a pdb, move it to another directory and return the filepath

    :param indir: location of pdb.gz files
    :param outdir: output location
    :return: filepath of .pdb
    """
    ## unzip pdb and move to outdir
    try:
        with gzip.open('%s/IMGT-%s.pdb.gz' % (indir, ID), 'rb') as f_in:
            # Check if the file is empty
            if f_in.seek(0, whence=2) == 0:
                raise Exception('File is empty')

        with gzip.open('%s/IMGT-%s.pdb.gz' % (indir, ID), 'rb') as f_in:
            with open('%s/%s.pdb' % (outdir, ID), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    except FileNotFoundError:
        print('ERROR TYPE #1: File not found. %s' % ID)

    return '%s/%s.pdb' % (outdir, ID)

def find_chains_MHCII(ID, pdb_file,pept_dist = 10):
    ''' Find the Alpha, Beta and Peptide chain in a pdb file. Also checks if there are heteroatoms in the binding
    cleft and checks for the proper size of the peptide chain.

    :param ID: (string) pdb ID
    :param pdb_file: (string) path to pdb file
    :param pept_dist: (int) distance in angstrom between the peptide chain and putative MHC chains to select them.
    :return: [Alpha_chain, Beta_chain, Peptide_chain] and the pdb object
    '''
    parser = PDBParser(QUIET=True)  # Create a parser object, used to read pdb files
    pdb = parser.get_structure(ID, pdb_file)

    # Remove weird chains
    for i in pdb.get_chains():
        for model in pdb:
            for chain in model:
                if chain.id in [' ', '', '1', '2', '3', '4', '5', '6', '7', '8', 'S']:
                    model.detach_child(chain.id)

    try:
    # Find most likely peptide chain: first chain to be 7 < len(chain) < 25
        pept_chain = []
        for chain in pdb.get_chains():
            if len(chain) > 6 and len(chain) < 26: # Is this chain between 7 and 25?
                heteroatoms = False
                for res in chain:
                    if res.id[0] != " " and res.id[0] != 'W': # Check if a res in this chain is a heteroatom
                        heteroatoms = True
                        break
                if heteroatoms == False: # If all residues are oke, add this to the list of peptide chains
                    pept_chain.append(chain.id)
        pept_chain = pept_chain[0] #Take the first pept chain. If there are multiple, they are probably duplicates

        if len(pept_chain) == 0:
            print('Could not find peptide chain')
            raise Exception
    except:
        print('Could not find peptide chain')
        raise Exception

    # Find peptide-structure contacts
    try:
        chain_cont = []
        for pept_res in pdb[0][pept_chain].get_residues(): #loop through all alpha carbons of the peptide residues
            for chain in pdb.get_chains(): #go through all chains and find residues of that chain that contact the peptide
                if chain.id != pept_chain: #dont measure distance between peptide-peptide residues of course
                    for res in chain:
                        if res.id[0] == ' ': #check if the residue is an amino acid and no heteroatom or water
                            try:
                                if res['CA'] - pept_res['CA'] < pept_dist: #distance between peptide and chain residue alpha carbon
                                    chain_cont.append(chain.id)
                            except:
                                pass
    except:
        print('Could not find peptide chain')
        raise Exception

    # Find the two chains that have the most contacts with the peptide. This should also filter out TCR chains
    if len(set(chain_cont)) >= 2:
        MHC_chains = sorted([ss for ss in set(chain_cont)], key=chain_cont.count, reverse=True)[:2]
        MHC_chains = sorted(MHC_chains) # sort alphabetically
        MHC_chains = MHC_chains + [pept_chain]
    else:
        print('Found >2 MHC chains')
        raise Exception

    for i in pdb.get_chains():
        for model in pdb:
            for chain in model:
                if chain.id not in MHC_chains:
                    model.detach_child(chain.id)

    return MHC_chains, pdb

def replace_chain_names(chains, pdb, replacement_chains=['M', 'N', 'P']):
    ''' Replace chain names by another chain name in a bio.pdb object

    :param chains: (list) chains to replace
    :param pdb: bio.pdb object
    :param replacement_chains: (list) replacement names (in order of chains to replace)
    :return: bio.pdb object with changed chain names
    '''

    for i in chains:
        for chain in pdb.get_chains():
            if chain.id == i:
                chain.id = replacement_chains[chains.index(i)]
    return pdb


def write_pdb(pdb, out_path, get_header_from=False):
    ''' Write bio.pdb object to file, can use the header of the original pdb (bio.pdb cant remember file headers)

    :param pdb: bio.pdb object
    :param out_path: (string)
    :param get_header_from: (string) get the header from another pdb file
    '''

    def get_head_and_remarks(pdb_file):
        ''' Get the head and remarks of an IMGT pdb file

        :param pdb_file: path to pdb file
        :return: (list) list of lines
        '''
        # Count until where the header and remarks last
        x = 0
        last_header_line = 0
        with open(pdb_file) as infile:
            header = []
            for line in infile:
                x += 1
                header.append(line[:-1])
                if line.startswith('REMARK'):
                    if last_header_line == 0:
                        last_header_line = x
                    last_header_line += 1
        header = [x for x in header[:last_header_line - 1] if
                  x != []]  # remove everything after the remarks and whitelines
        return header

    def line_prepender(filename, line):
        ''' Add a line in front of a file

        :param filename: (string) filepath
        :param line: (string) line to prepend
        '''
        with open(filename, 'r+') as f:
            content = f.read()
            f.seek(0, 0)
            f.write(line.rstrip('\r\n') + '\n' + content)

    # If the original pdb file path is given, use that header and paste it before the ATOM lines
    if get_header_from:
        header = get_head_and_remarks(get_header_from)  # get the header

    # Write pdb
    io = PDBIO()
    io.set_structure(pdb)
    io.save(out_path)

    if get_header_from: # Write the original header to pdb file
        for line in header[::-1]:
            line_prepender(out_path, line)

def parse_pMHCII_pdbs(ids_list):
    ''' Clean all MHCII pdb files that have been downloaded from IMGT

    :param ids_list: (list) list of MHCII PDB IDs. core.Database.IDs_list_MHCII
    :return: Writes all cleaned PDBs to the /PDBs/pMHCII/ dir
    '''
    # keep track of the failed PDBs
    bad_ids = []
    # set paths for in and out directories
    indir = PANDORA.PANDORA_data + '/PDBs/IMGT_retrieved/IMGT3DFlatFiles'
    outdir = PANDORA.PANDORA_data + '/PDBs/pMHCII'
    bad_dir = PANDORA.PANDORA_data + '/PDBs/Bad/pMHCII'

    for ID in ids_list:
        # ID = '1D5M'
        print('Parsing %s' %ID)
        try:
            # Unzip file (also check if the file is not empty) and save the path of this file
            pdb_file = unzip_pdb(ID, indir, outdir)
            # Find out which chains are the Alpha, Beta and Peptide chain
            MHC_chains, pdb = find_chains_MHCII(ID, pdb_file)
            # Rename the Alpha, Beta and Peptide chain to M,N,P respectively
            pdb = replace_chain_names(MHC_chains, pdb, ['M', 'N', 'P'])
            # Finally, write the cleaned pdb to the output dir. Keep the header of the original file.
            write_pdb(pdb, '%s/%s.pdb' %(outdir, ID), pdb_file)

        except: # If something goes wrong, append the ID to the bad_ids list
            bad_ids.append(ID)
            os.system('mv %s/%s.pdb %s/%s.pdb' %(outdir.replace(' ', '\\ '), ID, bad_dir.replace(' ', '\\ '), ID))
            pass

    return bad_ids


