import PANDORA
import pickle
from PANDORA.PMHC import PMHC
from PANDORA.Database import Database_functions
import os
import subprocess
from joblib import Parallel, delayed

class Database:

    def __init__(self):
        self.MHCI_data = {}
        self.MHCII_data = {}
        self.ref_MHCI_sequences = {}
        self.__IDs_list_MHCI = []
        self.__IDs_list_MHCII = []

    def download_data(self, data_dir = PANDORA.PANDORA_data, download = True):
        """download_data(self, data_dir = PANDORA.PANDORA_data, download = True)
        Download all MHC structures and get a two lists that contains all MHCI and MHCII IDs respectively"""

        if download:
            print('Downloading structures ...')
            Database_functions.download_unzip_imgt_structures(data_dir = data_dir, del_inn_files=True, del_kabat_files=True)
        self.__IDs_list_MHCI = Database_functions.download_ids_imgt('MH1', data_dir = data_dir, out_tsv='all_MHI_IDs.tsv')
        self.__IDs_list_MHCII = Database_functions.download_ids_imgt('MH2', data_dir = data_dir, out_tsv='all_MHII_IDs.tsv')


    def clean_MHCI_file(self, pdb_id, data_dir, remove_biopython_object):
        """ Clean all MHCI structures"""
        try: 
            templ = Database_functions.parse_pMHCI_pdb(pdb_id,
                                                    indir = data_dir + '/PDBs/IMGT_retrieved/IMGT3DFlatFiles',
                                                    outdir = data_dir + '/PDBs/pMHCI',
                                                    bad_dir = data_dir + '/PDBs/Bad/pMHCI',
                                                    remove_biopython_object=remove_biopython_object)
            if templ != None:
                #self.MHCI_data[pdb_id] = templ
                return (pdb_id, templ)

        except Exception as e:
            print('something went wrong parsing %s:' %pdb_id)
            print(e)

    def clean_MHCII_file(self, pdb_id, data_dir, remove_biopython_object):
        """ Clean all MHCII structures. Returns a list of bad PDBs"""
        try: 
            templ = Database_functions.parse_pMHCII_pdb(pdb_id,
                                                   indir = data_dir + '/PDBs/IMGT_retrieved/IMGT3DFlatFiles',
                                                   outdir = data_dir + '/PDBs/pMHCII',
                                                   bad_dir = data_dir + '/PDBs/Bad/pMHCII',
                                                   remove_biopython_object=remove_biopython_object)
            if templ != None:
                #self.MHCI_data[pdb_id] = templ
                return (pdb_id, templ)

        except Exception as e:
            print('something went wrong parsing %s' %pdb_id)
            print(e)

    def update_ref_sequences(self):
        """Downloads and parse HLA and other MHC sequences to compile reference fastas.
        Returns a dictionary that can be used to select the desired reference sequence"""
        self.ref_MHCI_sequences = Database_functions.generate_mhcseq_database()

    def construct_database(self, save=PANDORA.PANDORA_data + '/PANDORA_database.pkl', data_dir = PANDORA.PANDORA_data,
                           MHCI=True, MHCII=True, download=True,
                           update_ref_sequences=True, 
                           remove_biopython_objects = True,
                           n_jobs = -1):
        '''construct_database(self, save, data_dir = PANDORA.PANDORA_data, MHCI=True, MHCII=True, download=True, update_ref_sequences=True)
        Construct the database. Download, clean and add all structures

        Args:
            save (str/bool): Filename of database pkl object. If False, does not save the database pkl. 
                If a path is provided, saved the database .pkl to that path. Defaults to PANDORA.PANDORA_data + '/default/PANDORA_database.pkl'.
            data_dir (str): Path of data directory. Defaults to PANDORA.PANDORA_data.
            MHCI (bool): Parse data for MHCI. Defaults to True.
            MHCII (bool): Parse data for MHCII. Defaults to True.
            download (bool): If True, download the structures data from IMGT. Defaults to True.
            update_ref_sequences (bool): If True, downloads and parse reference sequence strcutres. Defaults to True
            remove_biopython_objects (bool): If True, removes the biopython pdb 
                objects from the template objects to make the database considerably lighter.
                Switch to False only if the biopython objects are necessary. Defaults to True.
            
        Returns: Database object

        '''
        # Download the data
        self.download_data(download = download, data_dir = data_dir)

        # Construct the MHCI database
        if MHCI:
            # Parse all MHCI files
            templates = Parallel(n_jobs = n_jobs)(delayed(self.clean_MHCI_file)(id, data_dir, remove_biopython_objects) for id in self.__IDs_list_MHCI)
            templates = [x for x in templates if x != None]
            self.MHCI_data = {key: value for (key, value) in templates}

        # Construct the MHCII database
        if MHCII:
            # Parse all MHCII files
            templates = Parallel(n_jobs = n_jobs)(delayed(self.clean_MHCII_file)(id, data_dir, remove_biopython_objects) for id in self.__IDs_list_MHCII)
            templates = [x for x in templates if x != None]
            self.MHCII_data = {key: value for key, value in templates}

        #Download and parse HLA and MHC sequences reference data
        if update_ref_sequences:
            self.update_ref_sequences()

        self.construct_both_blast_db()

        if save:
            self.save(save)

        print('Database correctly generated')

    def add_structure(self, id, allele_type, peptide = '', 
                      MHC_class = 'I', chain_seq = [], anchors = [], 
                      pdb_path = False, pdb = False, remove_biopython_object=True):
        ''' Add a single structure to the database

        Args:
            id: (str) PDB identifier
            allele_type: (lst) list of MHC alleles (or allele)
            peptide: (str) peptide sequence
            MHC_class: (str) either 'I' or 'II' denoting MHC class I and MHC class II respectively
            chain_seq: (lst) list of chain sequence(s) for the M and N (Alpha and Beta) chain respectively
            anchors: (lst) list of integers specifying which residue(s) of the peptide should be fixed as an anchor
                        during the modelling. MHC class I typically has 2 anchors, while MHC class II typically has 4.
            pdb_path: (str) path to pdb file
            pdb: (Bio.PDB) Biopython PBD object

        '''

        if not id:
            if not pdb_path or not pdb:
                raise ValueError('Structure id or path of .pdb files was not given. Enter value for id and pdb_path')
        # Add to MHCI data
        if MHC_class == 'I':
            self.MHCI_data[id] = PMHC.Template(id, allele_type, peptide, 
                                               MHC_class, chain_seq, anchors, 
                                               pdb_path, pdb, remove_biopython_object)
        # Add to MHCII data
        if MHC_class == 'II':
            self.MHCII_data[id] = PMHC.Template(id, allele_type, peptide, 
                                                MHC_class, chain_seq, anchors, 
                                                pdb_path, pdb, remove_biopython_object)


    def write_db_into_fasta(self, outfile):
        """
        Writes structure db into a fasta file (to be later used to build a blast database)

        Args:
            outfile (str): output file path.

        Returns:
            None.

        """

        sequences = []
        for template in self.MHCI_data.values():
            #Get Header and sequence
            header, seq = Database_functions.get_sequence_for_fasta(
                                        template, MHC_class='I', chain='M')
            #Keep only the G-domain
            seq = seq[PANDORA.MHCI_G_domain[0][0]:PANDORA.MHCI_G_domain[0][1]]
            #Append to the list
            sequences.append((header, seq))

        for template in self.MHCII_data.values():
            #Get Header and sequence
            header, seq = Database_functions.get_sequence_for_fasta(
                                        template, MHC_class='II', chain='M')
            #Keep only the G-domain
            seq = seq[PANDORA.MHCII_G_domain[0][0]:PANDORA.MHCII_G_domain[0][1]]
            #Append to the list
            sequences.append((header, seq))

            #Get Header and sequence
            header, seq = Database_functions.get_sequence_for_fasta(
                                        template, MHC_class='II', chain='N')
            #Keep only the G-domain
            seq = seq[PANDORA.MHCII_G_domain[1][0]:PANDORA.MHCII_G_domain[1][1]]
            #Append to the list
            sequences.append((header, seq))

        self.all_sequences = sequences

        with open(outfile, 'w') as outfasta:
            for sequence in sequences:
                header = sequence[0]
                seq = sequence[1]
                outfasta.write('> %s\n' %header)
                outfasta.write('\n'.join(seq[j:j+60] for j in range(0, len(seq), 60)) + '\n')



    def construct_blast_db(self, infile, outpath, db_name):
        """
        Construc blast database for seq based template selection

        Args:
            outpath (str, optional): Data dir folder. Defaults to PANDORA.PANDORA_data.
            db_name (str, optional): Name of the db folder and fasta file. Defaults to 'MHC_blast_db'.
        Returns:
            None.

        """
        # if not os.path.isdir(outpath):
        #     subprocess.check_call('mkdir %s' %outpath, shell=True)

        # out_fasta = outpath+'/'+db_name+'.fasta'
        # self.write_db_into_fasta(outfile=out_fasta)

        subprocess.check_call((' ').join(['makeblastdb','-dbtype','prot',
                                          '-in', infile,'-out',
                                          outpath + '/' + db_name]), shell=True)

    def construct_both_blast_db(self, data_dir=PANDORA.PANDORA_data):

        #Define db name and path
        db_name = 'templates_blast_db'
        outpath = data_dir + '/' + db_name
        out_fasta = outpath + '/'+ db_name +'.fasta'

        #Create db directory
        if not os.path.isdir(outpath):
            subprocess.check_call('mkdir %s' %outpath, shell=True)

        #Create .fasta for the db
        self.write_db_into_fasta(outfile=out_fasta)

        #Construct blast database for blast-based sequence-based template selection
        self.construct_blast_db(infile = out_fasta,
                                outpath=outpath,
                                db_name=db_name)

        #Define db name and path
        db_name = 'refseq_blast_db'
        outpath = data_dir + '/' + db_name
        out_fasta = outpath + '/' + db_name + '.fasta'

        #Create db directory
        if not os.path.isdir(outpath):
            subprocess.check_call('mkdir %s' %outpath, shell=True)

        #Create .fasta for the db
        command='cat %s/mhcseqs/Human_MHC_data.fasta %s/mhcseqs/NonHuman_MHC_data.fasta > %s' %(data_dir,
                                                                              data_dir,
                                                                              out_fasta)
        subprocess.check_call(command, shell=True)

        #Construct blast database for retriving mhc allele
        self.construct_blast_db(infile=out_fasta,
                                outpath=outpath,
                                db_name=db_name)



    def remove_structure(self, id =''):
        ''' Removes a structure (by id) from the database

        Args:
            id: (str) PDB ID

        '''

        # Remove structure from database
        self.MHCI_data.pop(id, None)
        self.MHCII_data.pop(id, None)

    def save(self, fn = PANDORA.PANDORA_data + '/PANDORA_database.pkl'):
        """Save the database as a pickle file

        :param fn: (str) pathname of file
        """
        with open(fn, "wb") as pkl_file:
            pickle.dump(self, pkl_file)

def load(file_name = PANDORA.PANDORA_data + '/PANDORA_database.pkl'):
    """Loads a pre-generated database


    Args:
        file_name (str): Dabase file name/path. 
            Defaults to PANDORA.PANDORA_data + '/PANDORA_database.pkl'.

    Returns:
        Database.Database: Database object.

    Example:
        >>> db = Database.load('MyDatabase.pkl')

    """
    try:
        with open(file_name, 'rb') as inpkl:
            db = pickle.load(inpkl)
        return db
    except FileNotFoundError:
        raise Exception('Database file not found. Are you sure you have it? If not, run Database.construct_database()')
