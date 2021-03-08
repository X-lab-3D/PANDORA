
import os
import pickle

import PANDORA
from PANDORA.PMHC import PMHC
from PANDORA.Database import Download_data
from PANDORA.Database import Parse_pMHC

class Database:
    #todo Integrate Anchor calculation with Rafaellas code to do all the anchor calcs while initiating the db
    def __init__(self):
        self.MHCI_data = {}
        self.MHCII_data = {}
        self.__IDs_list_MHCI = []
        self.__IDs_list_MHCII = []

    def download_data(self):
        """ Download all MHC structures and get a two lists that contains all MHCI and MHCII IDs respectively"""
        print('Downloading structures ...')
        # Download_data.download_unzip_imgt_structures(del_inn_files=True, del_kabat_files=True)
        self.__IDs_list_MHCI = Download_data.download_ids_imgt('MH1', out_tsv='all_MHI_IDs.tsv')
        self.__IDs_list_MHCII = Download_data.download_ids_imgt('MH2', out_tsv='all_MHII_IDs.tsv')


    def __clean_MHCI_files(self):
        """ Clean all MHCI structures"""
        return Parse_pMHC.parse_pMHCI_pdbs(self.__IDs_list_MHCI)

    def __clean_MHCII_files(self):
        """ Clean all MHCII structures. Returns a list of bad PDBs"""
        return Parse_pMHC.parse_pMHCII_pdbs(self.__IDs_list_MHCII)

    def construct_database(self, MHCI = True, MHCII = True):
        """ Construct a database

        :param MHCI: (bool) Calculate metadata for MHCI
        :param MHCII: (bool) Calculate metadata for MHCII
        :return: (dict) {id:MHC_structure}
        """

        # Download the data
        self.download_data()

        # Construct the MHCI database
        if MHCI:

            # Parse all MHCI files
            # self.__clean_MHCI_files()

            for id in self.__IDs_list_MHCI:
                try:
                    file = PANDORA.PANDORA_data + '/PDBs/pMHCI/' + id + '.pdb'
                    # Check if this file really exists
                    if os.path.isfile(file):
                        print('Adding %s' % id)
                        # Find the alleles
                        al = Parse_pMHC.get_chainid_alleles_MHCI(file)
                        alpha = [[k for k,v in i.items()] for i in [al['A'][i] for i in [i for i in al['A'].keys()]]]
                        a_allele = sum(alpha, [])
                        # Create MHC_structure object
                        self.MHCI_data[id] = PMHC.Template(id, allele=a_allele, pdb_path=file)
                except:
                    pass

        # Construct the MHCII database
        if MHCII:

            # Parse all MHCII files
            # self.__clean_MHCII_files()

            for id in self.__IDs_list_MHCII:
                try:
                    file = PANDORA.PANDORA_data + '/PDBs/pMHCII/' + id + '.pdb'
                    # Check if this file really exists
                    if os.path.isfile(file):
                        print('Adding %s' %id)
                        # Find the alleles
                        al = Parse_pMHC.get_chainid_alleles_MHCII(file)
                        alpha = sum([al['Alpha'][i] for i in [i for i in al['Alpha'].keys()]], [])
                        beta = sum([al['Beta'][i] for i in [i for i in al['Beta'].keys()]], [])
                        a_allele = list(set([alpha[i - 1] for i in range(3, int(len(alpha)), 4)]))
                        b_allele = list(set([beta[i - 1] for i in range(3, int(len(beta)), 4)]))

                        # Create MHC_structure object
                        self.MHCII_data[id] = PMHC.Template(id, allele=a_allele + b_allele, MHC_class= 'II', pdb_path=file)
                except:
                    pass


    def add_structure(self, PDB_id, allele, peptide = '', MHC_class = 'I', chain_seq = [], anchors = [], pdb_path = False, pdb = False):
        ''' Add a single structure to the database. Needs id and pdb_file as input. More can be given.
            Input is the same as the input from MHC_structure.template(). Allele needs to be given as input as well and
            is not fetched on the go. If no MHC class if specified, it will assume the structure is MHCI

        :param PDB_id: (string) PDB identifier
        :param allele: (list) list of MHC alleles (or allele)
        :param peptide: (string) peptide sequence
        :param MHC_class: (string) either 'I' or 'II' denoting MHC class I and MHC class II respectively
        :param chain_seq: (list) list of chain sequence(s) for the M and N (Alpha and Beta) chain respectively
        :param anchors: (list) list of integers specifying which residue(s) of the peptide should be fixed as an anchor
                        during the modelling. MHC class I typically has 2 anchors, while MHC class II typically has 4.
        :param pdb_path: (string) path to pdb file
        :param pdb: (Bio.PDB) Biopython PBD object
        '''

        if not PDB_id:
            if not pdb_path or not pdb:
                raise ValueError('Structure id or path of .pdb files was not given. Enter value for PDB_id and pdb_path')
        # Add to MHCI data
        if MHC_class == 'I':
            self.MHCI_data[id] = PMHC.Template(PDB_id, allele, peptide, MHC_class, chain_seq, anchors, pdb_path, pdb)
        # Add to MHCII data
        if MHC_class == 'II':
            self.MHCII_data[id] = PMHC.Template(PDB_id, allele, peptide, MHC_class, chain_seq, anchors, pdb_path, pdb)

    def remove_structure(self, id =''):
        ''' Removes a structure (by id) from the database

        :param id: (string) PDB ID
        '''
        # Remove structure from database
        self.MHCI_data.pop(id, None)
        self.MHCII_data.pop(id, None)

    # def save(self, fn = 'db.pkl'):
    #     """Save the database as a pickle file
    #
    #     :param fn: (string) pathname of file
    #     """
    #     pickle.dump(self, open(fn, "wb"))
    #
    # def load(cls, fn):
    #     return pickle.load( open(fn, "rb" ) )

    def save(self, fn = 'db.pkl'):
        """Save the database as a pickle file

        :param fn: (string) pathname of file
        """
        # pickle.dump(self, open(fn, "wb"))
        with open(fn, "wb") as dill_file:
            dill.dump(self, dill_file)

    def load(cls, fn):
        # return pickle.load( open(fn, "rb" ) )
        return dill.load(open(fn, 'rb'))

# db = Database()
# db.construct_database()
# db.save('Pandora_MHCI_and_MHCII_data.pkl')
# test_db = pickle.load( open("Pandora_MHCI_and_MHCII_data.pkl", "rb" ) )


