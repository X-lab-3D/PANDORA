
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1

import PANDORA
import os
import pickle

from PANDORA.tools import pdb_atom_contacts
from PANDORA.parsing import url_protocols
from PANDORA.parsing import structures_parser

class PMHC:

    def __init__(self, PDB_id, allele, peptide = '', MHC_class = 'I', chain_seq = [], anchors = []):
        ''' pMHC class. Acts as a parent class to Template and Target.

        :param PDB_id: (string) PDB identifier
        :param allele: (list) list of MHC alleles (or allele)
        :param peptide: (string) peptide sequence
        :param MHC_class: (string) either 'I' or 'II' denoting MHC class I and MHC class II respectively
        :param chain_seq: (list) list of chain sequence(s) for the M and N (Alpha and Beta) chain respectively
        :param anchors: (list) list of integers specifying which residue(s) of the peptide should be fixed as an anchor
                        during the modelling. MHC class I typically has 2 anchors, while MHC class II typically has 4.
        '''
        self.PDB_id = PDB_id
        self.MHC_class = MHC_class
        self.peptide = peptide
        self.chain_seq = chain_seq
        self.allele = allele
        self.anchors = anchors


class Template(PMHC):

    def __init__(self, PDB_id, allele, peptide = '', MHC_class = 'I', chain_seq = [], anchors = [], pdb_path = False, pdb = False):
        ''' Template structure class. This class holds all information of a template structure that is used for
            homology modelling. This class needs a PDB_id, allele and the path to a pdb file to work. (sequence info of
            the chains and peptide can be fetched from the pdb)

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
        super().__init__(PDB_id, allele, peptide, MHC_class, chain_seq, anchors)
        self.pdb_path = pdb_path
        self.pdb = pdb
        self.chain_contacts = False

        if not pdb_path or not pdb: # If the path to a pdb file or a Bio.PDB object is given, parse the pdb
            self.parse_pdb()

    def parse_pdb(self):
        '''Loads pdb from path, updates self.pdb field and self.chain_seq/self.peptide if they were empty'''

        if self.pdb_path and not self.pdb: #if there is a path to a pdb provided and there is not already a self.pdb...
            parser = PDBParser(QUIET=True)  # Create a parser object, used to read pdb files
            self.pdb = parser.get_structure('MHC', self.pdb_path) #Create Bio.PDB object

        # If the chains or peptide are not given by the user, fetch them from the pdb
        # Get the chain sequences
        chain_seqs = [seq1(''.join([res.resname for res in chain])) for chain in self.pdb.get_chains()]
        # Update chain and peptide fields if emtpy
        if self.MHC_class == 'I':
            if not self.chain_seq:
                self.chain_seq = [chain_seqs[0]]
            if not self.peptide:
                self.peptide = chain_seqs[-1]
        if self.MHC_class == 'II':
            if not self.chain_seq:
                self.chain_seq = chain_seqs[:2]
            if not self.peptide:
                self.peptide = chain_seqs[-1]

    def info(self):
        """ Print the basic info of this structure

        """
        print('This is a %s structure.' %(type(self).__name__))
        print('ID: %s' %self.PDB_id)
        print('Type: MHC class %s' %self.MHC_class)
        print('Alleles: %s' % self.allele)
        if len(self.chain_seq) > 0:
            print('Alpha chain length: %s' %len(self.chain_seq[0]))
        if len(self.chain_seq) > 1:
            print('Beta chain length: %s' %len(self.chain_seq[1]))
        print('Peptide length: %s' %len(self.peptide))
        if len(self.chain_seq) > 0:
            print('Alpha chain: %s' % self.chain_seq[0])
        if len(self.chain_seq) > 1:
            print('Beta chain: %s' % self.chain_seq[1])
        print('Peptide: %s' % self.peptide)
        print('Anchors: %s' %self.anchors)
        if self.pdb_path:
            print('Path to PDB file: %s' %self.pdb_path)
        if not self.pdb:
            print('PDB structure: no PDB structure provided')
        else:
            print('PDB structure:')
            for (k, v) in self.pdb.header.items():
                print('\t'+k + ':', v)


    def calc_contacts(self):
        if self.pdb:
            self.chain_contacts = pdb_atom_contacts.Contacts(self.pdb)
        else:
            raise Exception('Provide a PDB structure to the Template object first')

    def calc_anchors(self):
        pass


class Target(PMHC):

    def __init__(self, PDB_id, peptide, allele, MHC_class = 'I', chain_seq = [], anchors = [], use_template = False):
        ''' Target structure class. This class needs an ID (preferably a PDB ID), allele and pepide information.

        :param PDB_id: (string) PDB identifier
        :param allele: (list) list of MHC alleles (or allele)
        :param peptide: (string) peptide sequence
        :param MHC_class: (string) either 'I' or 'II' denoting MHC class I and MHC class II respectively
        :param chain_seq: (list) list of chain sequence(s) for the M and N (Alpha and Beta) chain respectively
        :param anchors: (list) list of integers specifying which residue(s) of the peptide should be fixed as an anchor
                        during the modelling. MHC class I typically has 2 anchors, while MHC class II typically has 4.
        :param use_template: (string) or (Bio.PDB) The user can specify that PANDORA uses a certain structure as
                        template. This can be provided as the path to a pdb file or a Bio.PDB object. The pdb must
                        contain of a Alpha (and Beta for MHC Type II) chain and peptide, without extra chains like TCR
                        or small molecules.

        '''
        super().__init__(PDB_id, peptide, allele, MHC_class, chain_seq, anchors)
        self.use_template = use_template
        self.initial_model = False
        self.anchor_contacts = False

    def info(self):
        """ Print the basic info of this structure

        """
        print('This is a %s structure.' % (type(self).__name__))
        print('ID: %s' % self.PDB_id)
        print('Type: MHC class %s' % self.MHC_class)
        print('Alleles: %s' % self.allele)
        if len(self.chain_seq) > 0:
            print('Alpha chain length: %s' % len(self.chain_seq[0]))
        if len(self.chain_seq) > 1:
            print('Beta chain length: %s' % len(self.chain_seq[1]))
        print('Peptide length: %s' % len(self.peptide))
        if len(self.chain_seq) > 0:
            print('Alpha chain: %s' % self.chain_seq[0])
        if len(self.chain_seq) > 1:
            print('Beta chain: %s' % self.chain_seq[1])
        print('Peptide: %s' % self.peptide)
        print('Anchors: %s' % self.anchors)
        if self.use_template:
            print('Using template %s for homology modelling' %self.use_template) #todo actually build this functionality
        if self.initial_model:
            print('An initial model has been provided.')


    def calc_anchor_contacts(self):
        if self.initial_model and self.anchors:
            self.anchor_contacts = pdb_atom_contacts.Contacts(self.pdb, anchors=self.anchors)
        else:
            raise Exception('Provide an initial model (.ini PDB) and anchor positions to the Target object first')

    def calc_anchors(self):
        pass

# test = Template('XXXX', allele = ['allele1','allele2'], MHC_class= 'II', pdb_path =  PANDORA.PANDORA_data + '/PDBs/pMHCII/1IAK.pdb')
#
# tar = Target('XXXX', 'NNNNN', [])


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
        # url_protocols.download_unzip_imgt_structures(del_inn_files=True, del_kabat_files=True)
        self.__IDs_list_MHCI = url_protocols.download_ids_imgt('MH1', out_tsv='all_MH1_IDs.tsv')
        self.__IDs_list_MHCII = url_protocols.download_ids_imgt('MH2', out_tsv='all_MH1I_IDs.tsv')


    def __clean_MHCI_files(self):
        """ Clean all MHCI structures"""
        structures_parser.parse_pMHCI_pdbs(self.__IDs_list_MHCI)

    def __clean_MHCII_files(self):
        """ Clean all MHCII structures"""
        structures_parser.parse_pMHCII_pdbs(self.__IDs_list_MHCII)

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
            self.__clean_MHCI_files()

            for id in self.__IDs_list_MHCI:
                file = PANDORA.PANDORA_data + '/PDBs/pMHCI/' + id + '.pdb'
                # Check if this file really exists
                if os.path.isfile(file):
                    # Find the alleles
                    al = structures_parser.get_chainid_alleles_MHCI(file)
                    alpha = sum([al['Alpha'][i] for i in [i for i in al['Alpha'].keys()]], [])
                    a_allele = list(set([alpha[i - 1] for i in range(3, int(len(alpha)), 4)]))
                    # Create MHC_structure object
                    self.MHCI_data[id] = Template(id, allele=a_allele, pdb_path=file)

        # Construct the MHCII database
        if MHCII:

            # Parse all MHCII files
            self.__clean_MHCII_files()

            for id in self.__IDs_list_MHCII:
                file = PANDORA.PANDORA_data + '/PDBs/pMHCII/' + id + '.pdb'
                # Check if this file really exists
                if os.path.isfile(file):
                    # Find the alleles
                    al = structures_parser.get_chainid_alleles_MHCII(file)
                    alpha = sum([al['Alpha'][i] for i in [i for i in al['Alpha'].keys()]], [])
                    beta = sum([al['Beta'][i] for i in [i for i in al['Beta'].keys()]], [])
                    a_allele = list(set([alpha[i - 1] for i in range(3, int(len(alpha)), 4)]))
                    b_allele = list(set([beta[i - 1] for i in range(3, int(len(beta)), 4)]))
                    # Create MHC_structure object
                    self.MHCII_data[id] = Template(id, allele=a_allele + b_allele, MHC_class= 'II', pdb_path=file)


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
            self.MHCI_data[id] = Template(PDB_id, allele, peptide, MHC_class, chain_seq, anchors, pdb_path, pdb)
        # Add to MHCII data
        if MHC_class == 'II':
            self.MHCII_data[id] = Template(PDB_id, allele, peptide, MHC_class, chain_seq, anchors, pdb_path, pdb)

    def remove_structure(self, id =''):
        ''' Removes a structure (by id) from the database

        :param id: (string) PDB ID
        '''
        # Remove structure from database
        self.MHCI_data.pop(id, None)
        self.MHCII_data.pop(id, None)

    def save(self, fn):
        """Save the database as a pickle file

        :param fn: (string) pathname of file
        """
        db = open(fn, "wb")
        pickle.dump(self, db)
        db.close()

    def load(fn):
        """Load the database from a saved pickle file

        :param fn: (string) pathname of file
        :return: local_database object
        """
        pkl = open(fn, 'rb')
        db = pickle.load(pkl)
        pkl.close()
        return db


# db = Database()
# # db.download_data()
#
# db.construct_database(MHCI=False)
#
# struc = db.MHCII_data['1A6A']
# struc.PDB_id


class Pandora:

    def __init__(self):
        # self.target
        # self.template
        # self.database
        # self.alignment
        # self.ini_modeller_script
        # self.initial_model
        # self.modeller_script
        # self.results
        pass

    def find_template(self):
        pass

    def find_anchors(self):
        pass

    def align(self):
        pass

    def write_ini_script(self):
        pass

    def run_modeller(self):
        pass

    def anchor_contacts(self):
        pass

    def write_modeller_script(self):
        pass

    def model(self):
        pass











