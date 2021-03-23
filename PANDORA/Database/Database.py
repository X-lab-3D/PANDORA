import PANDORA
import dill
from PANDORA.PMHC import PMHC
from PANDORA.Database import Database_functions

class Database:
    #todo Integrate Anchor calculation with Rafaellas code to do all the anchor calcs while initiating the db
    def __init__(self):
        self.MHCI_data = {}
        self.MHCII_data = {}
        self.__IDs_list_MHCI = []
        self.__IDs_list_MHCII = []

    def download_data(self, data_dir = PANDORA.PANDORA_data, download = True):
        """ Download all MHC structures and get a two lists that contains all MHCI and MHCII IDs respectively"""
        print('Downloading structures ...')
        if download:
            Database_functions.download_unzip_imgt_structures(data_dir = data_dir, del_inn_files=True, del_kabat_files=True)
        self.__IDs_list_MHCI = Database_functions.download_ids_imgt('MH1', data_dir = data_dir, out_tsv='all_MHI_IDs.tsv')
        self.__IDs_list_MHCII = Database_functions.download_ids_imgt('MH2', data_dir = data_dir, out_tsv='all_MHII_IDs.tsv')


    def __clean_MHCI_file(self, pdb_id, data_dir):
        """ Clean all MHCI structures"""
        return Database_functions.parse_pMHCI_pdb(pdb_id,
                                                   indir = data_dir + '/PDBs/IMGT_retrieved/IMGT3DFlatFiles',
                                                   outdir = data_dir + '/PDBs/pMHCI',
                                                   bad_dir = data_dir + '/PDBs/Bad/pMHCI')

    def __clean_MHCII_file(self, pdb_id, data_dir):
        """ Clean all MHCII structures. Returns a list of bad PDBs"""
        return Database_functions.parse_pMHCII_pdb(pdb_id,
                                                   indir = data_dir + '/PDBs/IMGT_retrieved/IMGT3DFlatFiles',
                                                   outdir = data_dir + '/PDBs/pMHCII',
                                                   bad_dir = data_dir + '/PDBs/Bad/pMHCII')

    def construct_database(self, save, data_dir = PANDORA.PANDORA_data, MHCI=True, MHCII=True, download=True):
        ''' Construct the database. Download, clean and add all structures

        Args:
            save: (string/bool) Filename of database, can be False if you don't want to save the database
            data_dir: (string) Path of data directory
            MHCI: (bool) Parse data for MHCI
            MHCII: (bool) Parse data for MHCII
            download: (bool) Download the data? If you already downloaded the data, this can be set to False

        Returns: Database object

        '''
        # Download the data
        self.download_data(download = download, data_dir = data_dir)

        # Construct the MHCII database
        if MHCI:
            # Parse all MHCI files
            for id in self.__IDs_list_MHCI:
                try:
                    templ = self.__clean_MHCI_file(pdb_id = id, data_dir = data_dir)
                    if templ != None:
                        self.MHCI_data[id] = templ
                except:
                    pass

        # Construct the MHCII database
        if MHCII:
            # Parse all MHCII files
            for id in self.__IDs_list_MHCII:
                try:
                    templ = self.__clean_MHCII_file(pdb_id = id, data_dir = data_dir)
                    if templ != None:
                        self.MHCII_data[id] = templ
                except:
                    pass

        if save:
            self.save(save)

    def add_structure(self, id, allele_type, peptide = '', MHC_class = 'I', chain_seq = [], anchors = [], pdb_path = False, pdb = False):
        ''' Add a single structure to the database

        Args:
            id: (string) PDB identifier
            allele_type: (list) list of MHC alleles (or allele)
            peptide: (string) peptide sequence
            MHC_class: (string) either 'I' or 'II' denoting MHC class I and MHC class II respectively
            chain_seq: (list) list of chain sequence(s) for the M and N (Alpha and Beta) chain respectively
            anchors: (list) list of integers specifying which residue(s) of the peptide should be fixed as an anchor
                        during the modelling. MHC class I typically has 2 anchors, while MHC class II typically has 4.
            pdb_path: (string) path to pdb file
            pdb: (Bio.PDB) Biopython PBD object

        '''

        if not id:
            if not pdb_path or not pdb:
                raise ValueError('Structure id or path of .pdb files was not given. Enter value for id and pdb_path')
        # Add to MHCI data
        if MHC_class == 'I':
            self.MHCI_data[id] = PMHC.Template(id, allele_type, peptide, MHC_class, chain_seq, anchors, pdb_path, pdb)
        # Add to MHCII data
        if MHC_class == 'II':
            self.MHCII_data[id] = PMHC.Template(id, allele_type, peptide, MHC_class, chain_seq, anchors, pdb_path, pdb)

    def remove_structure(self, id =''):
        ''' Removes a structure (by id) from the database

        Args:
            id: (string) PDB ID

        '''

        # Remove structure from database
        self.MHCI_data.pop(id, None)
        self.MHCII_data.pop(id, None)

    def save(self, fn = 'db.pkl'):
        """Save the database as a pickle file

        :param fn: (string) pathname of file
        """
        with open(fn, "wb") as dill_file:
            dill.dump(self, dill_file)

    def load(cls, fn):
        return dill.load(open(fn, 'rb'))

#
