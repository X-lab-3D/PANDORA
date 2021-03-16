import os
import urllib.request
import urllib.parse
from string import ascii_uppercase
from copy import deepcopy
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import parse_pdb_header
import gzip
import shutil
import PANDORA
from PANDORA.Contacts import Contacts


def download_unzip_imgt_structures(del_inn_files = True, del_kabat_files = True):
    '''
    Downloads the complete structural dataset
    from IMGT database: http://www.imgt.org/download/3Dstructure-DB/IMGT3DFlatFiles.tgz
    The only two arguments delete every non-PDB structure file

    Args:
        del_inn_files(bool) : if True (default) deletes all inn files

        del_kabat_files(bool) : if True (default) deletes all kabat files

    '''

    # Changing working directory
    os.chdir(PANDORA.PANDORA_data + '/PDBs/IMGT_retrieved/')
    # Downloading IMGT dataset
    os.system('wget http://www.imgt.org/download/3Dstructure-DB/IMGT3DFlatFiles.tgz')
    # Uncompressing
    os.system('gunzip IMGT3DFlatFiles.tgz')
    os.system('tar -xvf IMGT3DFlatFiles.tar')

    try:
        os.system('rm IMGT3DFlatFiles.tgz')
    except:
        pass
    os.system('rm IMGT3DFlatFiles.tar')
    # Removing non-PDB files
    if del_inn_files:
        os.system('rm IMGT3DFlatFiles/*.inn.gz')
    if del_kabat_files:
        os.system('rm IMGT3DFlatFiles/*.prot.gz')
    #os.chdir('../../../../')

def download_ids_imgt(ReceptorType, out_tsv = False):
    '''
    Querys IMGT with the ReceptorType for PDBs.
    Returns the list of IDs provided by IMGT.

    Args:
        ReceptorType(str) : Receptor query for IMGT

        out_tsv(False or str) : if not False, produces a tsv file names as out_tsv

    Example:
    >>> pMHCI_list = download_ids_imgt('peptide/MH1', 'peptide_MHCI.tsv')

    '''

    #out_tsv = 'pMHCI_IDs_alleles_from_IMGT.tsv'
    params = { 'ReceptorType' : ReceptorType,
             'type-entry': 'PDB'}

    url = "http://www.imgt.org/3Dstructure-DB/cgi/3Dquery.cgi"

    data = urllib.parse.urlencode(params)
    data = data.encode('ascii') # data should be bytes
    req = urllib.request.Request(url, data)

    with urllib.request.urlopen(req) as response:
        text = response.read().decode('utf-8')
        text = text.splitlines()

    IDs_list = []
    IDs_list = [x for x in text if 'href' in x and 'pdbcode' in x]
    IDs_list = [x.split('"') for x in IDs_list]
    IDs_list = [x[3][-4:] for x in IDs_list]

    if out_tsv:
        outfile = open(PANDORA.PANDORA_data + '/csv_pkl_files/' + out_tsv, 'w')
        outfile.write(ReceptorType + ' IMGT IDs\n')
        for ID in IDs_list:
            outfile.write(ID + '\n')
        outfile.close()

    return IDs_list


def get_chainid_alleles_MHCI(pdbf):
    '''
    Takes as input an IMGT preprocessed PDB file of p:MHC I.
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

    ### Extracting MHC I Alpha chains
    mhc_a = {}  # MHC I Alpha
    for chain in chains:
        try:
            if chains[chain][1][3] == 'I-ALPHA':
                mhc_a[chain] = chains[chain]
        except:
            pass

    ### Extracting alleles info
    mhc_a_alleles = {}
    for chain in mhc_a:
        G_dom_alleles = {'G-ALPHA1': [], 'G-ALPHA2': []}
        key = False
        for row in mhc_a[chain]:
            if row[0] == 'G-DOMAIN':
                try:
                    if row[3] == 'description' and row[4] == 'G-ALPHA1':
                        key = 'G-ALPHA1'
                    elif row[3] == 'description' and row[4] == 'G-ALPHA2':
                        key = 'G-ALPHA2'
                    elif key:
                        if row[2] == 'gene' and row[3] == 'and' and row[4] == 'allele':
                            G_dom_alleles[key] += row[5:]
                        else:
                            key = False
                except IndexError:
                    pass
        mhc_a_alleles[chain] = deepcopy(G_dom_alleles)

    mhc_a_alleles_percs = {}
    for chain in mhc_a_alleles:
        mhc_a_alleles_percs[chain] = {}
        for key in mhc_a_alleles[chain]:
            mhc_a_alleles_percs[chain][key] = {}
            ### Allele info are always given with four elements: Gender, Spieces, Allele, Percentage
            for block in range(int(len(mhc_a_alleles[chain][key]) / 4)):
                allele = mhc_a_alleles[chain][key][2 + (4 * block)]
                perc = float(
                    mhc_a_alleles[chain][key][3 + (4 * block)].replace('(', '').replace('%)', '').replace(',', ''))
                mhc_a_alleles_percs[chain][key][allele] = perc

    return mhc_a_alleles_percs


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

    return {'Alpha': mhc_a_alleles, 'Beta': mhc_b_alleles}
#
# chains = MHC_chains, pdb, ['M', 'N', 'P']

def get_resolution(pdbf):
    ''' Returns the resolution in Angstrom from the given pdb

    Args:
        pdbf (str): path to the pdb file

    Returns:
        resolution (float): resolution of the model, in Angstrom

    '''
    header = parse_pdb_header(pdbf)
    resolution = header['resolution']
    return resolution



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
                # print(chain.id)
                chain.id = replacement_chains[chains.index(i)]

    return pdb

def renumber(pdb):
    ''' Renumbers the pdb. Each chain starts at 1

    :param pdb: Bio.PDb object
    :return: Bio.PDb object with renumbered residues
    '''

    for chain in pdb.get_chains():
        nr = 1
        for res in chain:
            res.id = ('X', nr, res.id[2])
            nr += 1
    for chain in pdb.get_chains():
        for res in chain:
            res.id = (' ', res.id[1], ' ')

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


def find_peptide_chain(pdb, min=6, max=26):
    ''' Find the pdb chain that is most likely the peptide based on its size

    Args:
        pdb: Bio.PDB object
        min: minimal peptide length to consider
        max: maximal peptide length to consider

    Returns: (string) Most likely chain that is the peptide

    '''
    # Find most likely peptide chain: first chain to be 7 < len(chain) < 25
    pept_chain = []
    for chain in pdb.get_chains():
        if len(chain) > min and len(chain) < max:  # Is this chain between 7 and 25?
            heteroatoms = False
            for res in chain:
                if res.id[0] != " " and res.id[0] != 'W':  # Check if a res in this chain is a heteroatom
                    heteroatoms = True
                    break
            if heteroatoms == False:  # If all residues are oke, add this to the list of peptide chains
                pept_chain.append(chain.id)
    pept_chain = pept_chain[0]  # Take the first pept chain. If there are multiple, they are probably duplicates

    if len(pept_chain) == 0:
        print('Could not find peptide chain')
        raise Exception

    return pept_chain

def remove_irregular_chains(pdb, chains_to_keep):
    ''' Removes all chains that you don't specify to keep

    Args:
        pdb: Bio.PDB object
        chains_to_keep: list of strings: ['A', 'C']

    Returns: Bio.PDB object

    '''
    for _ in range(len([c for c in pdb.get_chains()])):
        for i in pdb.get_chains():
            for model in pdb:
                for chain in model:
                    if chain.id not in chains_to_keep or chain.id in [' ','','  ','_'] + list(range(1,20)):
                        model.detach_child(chain.id)

    return pdb



def find_chains_MHCI(pdb, pept_chain):
    ''' Find the MHCI chains

    Args:
        pdb: Bio.PDB objet

    Returns: list of chains

    '''

    # Find contacts between peptide chain and other chains
    c = [i for i in Contacts.Contacts(pdb).chain_contacts if i[1] == pept_chain or i[5] == pept_chain]
    # Make a list of all chains that contact the peptide chain
    chain_cont = [i for i in sum([[i[1],i[5]] for i in c], []) if i != pept_chain]
    # Make sure the chain is longer than 120 residues. This prevents selecting e.g. two peptides in the binding groove
    chain_cont = [i for i in chain_cont if i in [c.id for c in pdb.get_chains() if len(c) > 120]]

    # Find the two chains that have the most contacts with the peptide. This should also filter out TCR chains
    if len(set(chain_cont)) >= 1:
        MHC_chains = sorted([ss for ss in set(chain_cont)], key=chain_cont.count, reverse=True)[0]
        MHC_chains = [MHC_chains, pept_chain]
    else:
        print('Found >1 MHC I chains')
        raise Exception

    return MHC_chains

def find_chains_MHCII(pdb, pept_chain):
    ''' Find the MHCI chains

    Args:
        pdb: Bio.PDB objet

    Returns: list of chains

    '''

    # Find contacts between peptide chain and other chains
    c = [i for i in Contacts.Contacts(pdb).chain_contacts if i[1] == pept_chain or i[5] == pept_chain]
    # Make a list of all chains that contact the peptide chain
    chain_cont = [i for i in sum([[i[1],i[5]] for i in c], []) if i != pept_chain]
    # Make sure the chain is longer than 120 residues. This prevents selecting e.g. two peptides in the binding groove
    chain_cont = [i for i in chain_cont if i in [c.id for c in pdb.get_chains() if len(c) > 120]]

    # Find the two chains that have the most contacts with the peptide. This should also filter out TCR chains
    if len(set(chain_cont)) >= 2:
        MHC_chains = sorted([ss for ss in set(chain_cont)], key=chain_cont.count, reverse=True)[:2]
        MHC_chains = sorted(MHC_chains)
        MHC_chains = MHC_chains + [pept_chain]
    else:
        print('Found >2 MHC I chains')
        raise Exception

    return MHC_chains


def check_pMHC(pdb):
    ''' Tests parsed pMHC structures: chain numbering, naming and length

    Args:
        pdb: Bio.PDB object

    Returns: Bool

    '''
    requirements = [False, False, False, True]

    chains = [i.id for i in pdb.get_chains()]
    chain_len = [len(i) for i in pdb.get_chains()]

    # 1. Check chain names and the number of chains
    if len(chains) == 2:
        if chains[0] == 'M' and chains[1] == 'P':
            requirements[0] = True
    elif len(chains) == 3:
        if chains[0] == 'M' and chains[1] == 'N' and chains[2] == 'P':
            requirements[0] = True

    # 2. Check M,N chain length
    if len(chains) == 2:
        if chain_len[0] > 120:
            requirements[1] = True
    elif len(chains) == 3:
        if chain_len[0] > 120 and chain_len[1] > 120:
            requirements[1] = True

    # 3. Check peptide length
    if chain_len[-1] > 6 and chain_len[-1] < 26:
        requirements[2] = True

    # 4. Check numbering. Does every chain start at 1 and end at its chain length?
    for c in pdb.get_chains():
        res = [r.id for r in c]
        if not res[0][1] == 1:
            requirements[3] = False
        if not res[-1][1] == len(c):
            requirements[3] = False

    if all(requirements):
        return True
    else:
        return False

def log(ID, error, logfile, verbose=True):
    ''' Keeps track of what goes wrong while parsing

    Args:
        ID: (string) PDB ID
        error: (string) error to append to log file
        logfile: (string) path to logfile
        verbose: (Bool) print error?
    '''

    # Create log file
    if not os.path.exists(logfile):
        with open(logfile, 'w') as f:
            f.write('ID,error\n')

    if verbose:
        print('\t' + error)
    with open(logfile, 'a') as f:
        f.write('%s,%s\n' % (ID, error))

# ID = IDs_list_MHCI[0]
def parse_pMHCI_pdbs(ids_list):
    ''' Clean all MHCI pdb files that have been downloaded from IMGT

    :param ids_list: (list) list of MHCI PDB IDs. core.Database.IDs_list_MHCI
    :return: Writes all cleaned PDBs to the /PDBs/pMHCI/ dir
    '''

    # set paths for in and out directories
    indir = PANDORA.PANDORA_data + '/PDBs/IMGT_retrieved/IMGT3DFlatFiles'
    outdir = PANDORA.PANDORA_data + '/PDBs/pMHCI'
    bad_dir = PANDORA.PANDORA_data + '/PDBs/Bad/pMHCI'
    logfile = os.path.dirname(bad_dir) + '/log_MHCI.csv'

    for ID in ids_list:
        print('Parsing %s' % ID)
        try:
            # Unzip file (also check if the file is not empty) and save the path of this file
            pdb_file = unzip_pdb(ID, indir, outdir)
            pdb = PDBParser(QUIET=True).get_structure('MHCI', pdb_file)
            # Remove waters and irregular chains
            pdb = remove_irregular_chains(pdb, list(ascii_uppercase))

            try:                # Find the peptide chain
                pept_chain = find_peptide_chain(pdb)
            except:
                log(ID, 'Could not find a suitable peptide chain', logfile)
                raise Exception

            try:                 # Find out which chains are the Alpha and Peptide chain
                MHC_chains = find_chains_MHCI(pdb, pept_chain)
            except:
                log(ID, 'Could not locate Alpha chain', logfile)
                raise Exception

            try:                 # Reformat chains
                pdb = remove_irregular_chains(pdb, MHC_chains)  # Remove all other chains from the PBD that we dont need
                pdb = replace_chain_names(MHC_chains, pdb,['M', 'P'])  # Rename chains to M,P
                pdb = renumber(pdb)  # Renumber from 1
            except:
                log(ID, 'Could not reformat structure', logfile)
                raise Exception

            if not check_pMHC(pdb):
                log(ID, 'Structure did not pass the test.', logfile)
                raise Exception

            # Finally, write the cleaned pdb to the output dir. Keep the header of the original file.
            write_pdb(pdb, '%s/%s.pdb' % (outdir, ID), pdb_file)

        except:  # If something goes wrong, append the ID to the bad_ids list
            os.system('mv %s/%s.pdb %s/%s.pdb' % (outdir, ID, bad_dir, ID))



def parse_pMHCII_pdbs(ids_list):
    ''' Clean all MHCII pdb files that have been downloaded from IMGT

    :param ids_list: (list) list of MHCI PDB IDs. core.Database.IDs_list_MHCII
    :return: Writes all cleaned PDBs to the /PDBs/pMHCII/ dir
    '''
    # set paths for in and out directories
    indir = PANDORA.PANDORA_data + '/PDBs/IMGT_retrieved/IMGT3DFlatFiles'
    outdir = PANDORA.PANDORA_data + '/PDBs/pMHCII'
    bad_dir = PANDORA.PANDORA_data + '/PDBs/Bad/pMHCII'
    logfile = os.path.dirname(bad_dir) + '/log_MHCII.csv'

    for ID in ids_list:
        print('Parsing %s' % ID)
        try:
            # Unzip file (also check if the file is not empty) and save the path of this file
            pdb_file = unzip_pdb(ID, indir, outdir)
            pdb = PDBParser(QUIET=True).get_structure('MHCII', pdb_file)
            # Remove waters and irregular chains
            pdb = remove_irregular_chains(pdb, list(ascii_uppercase))

            try:                # Find the peptide chain
                pept_chain = find_peptide_chain(pdb)
            except:
                log(ID, 'Could not find a suitable peptide chain', logfile)
                raise Exception

            try:                 # Find out which chains are the Alpha and Peptide chain
                MHC_chains = find_chains_MHCII(pdb, pept_chain)
            except:
                log(ID, 'Could not locate Alpha chain', logfile)
                raise Exception

            try:                 # Reformat chains
                pdb = remove_irregular_chains(pdb, MHC_chains)  # Remove all other chains from the PBD that we dont need
                pdb = replace_chain_names(MHC_chains, pdb, ['M', 'N', 'P'])   # Rename chains to M,N,P
                pdb = renumber(pdb)  # Renumber from 1
            except:
                log(ID, 'Could not reformat structure', logfile)
                raise Exception

            if not check_pMHC(pdb):
                log(ID, 'Structure did not pass the test.', logfile)
                raise Exception

            # Finally, write the cleaned pdb to the output dir. Keep the header of the original file.
            write_pdb(pdb, '%s/%s.pdb' % (outdir, ID), pdb_file)

        except:  # If something goes wrong, append the ID to the bad_ids list
            os.system('mv %s/%s.pdb %s/%s.pdb' % (outdir, ID, bad_dir, ID))