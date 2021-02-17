
from Bio.PDB import PDBParser
import PANDORA
from Bio.SeqUtils import seq1

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
        pass

    def __calc_contacts(self):
        pass

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
        pass

    def __calc_anchor_contacts(self):
        pass

    def calc_anchors(self):
        pass

test = Template('XXXX', allele = 'allele_1', MHC_class= 'II', pdb_path =  PANDORA.PANDORA_data + '/PDBs/pMHCII/1IAK.pdb')
























