from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO
import os
from Bio.PDB import PDBIO
from pdb2sql import StructureSimilarity

class Model:

    def __init__(self, target, output_dir, model_path='', pdb=False, molpdf=0, dope=0):
        '''

        Args:
            target: Target object
            output_dir: (string) output directory
            model_path: (string) path to hypothetical model
            pdb:  Bio.PDB object of the hypothetical model
            molpdf: (float) molpdf score
            dope:  (float) DOPE score
        '''


        self.target = target
        self.model_path = model_path
        self.moldpf = molpdf
        self.dope = dope
        self.output_dir = output_dir

        # Check if the user gave either the path to the model pdb or the pdb itself.
        if self.model_path == '' and not pdb:
            raise Exception('Provide the path to a model structure or a Bio.PDB object')
        # If there is a model path and no pdb, parse the pdb structure from that path.
        if not pdb:
            self.pdb = PDBParser(QUIET=True).get_structure(self.target.id, self.model_path)


    def calc_LRMSD(self, reference_pdb, atoms = ['C', 'CA', 'N', 'O']):
        ''' Calculate the L-RMSD between the decoy and reference structure (ground truth)

        Args:
            template_pdb: Bio.PDB object

        Returns: (float) L-RMSD
        '''

        # load target pdb
        if isinstance(reference_pdb, str):  # if its a string, it should be the path of the pdb, then load pdb first
            ref = PDBParser(QUIET=True).get_structure('MHC', reference_pdb)
        else:
            ref = template_pdb

        # pdb2sql needs 1 big chain and 1 ligand chain with correct numbering, for MHCII, this means merging the chains.
        homogenize_pdbs(self.pdb, ref, self.output_dir)

        # Calculate l-rmsd between decoy and reference with pdb2sql
        sim = StructureSimilarity('%s/decoy.pdb' % (self.output_dir.replace(' ', '\\ ')), '%s/ref.pdb' % (self.output_dir))
        # self.lrmsd = sim.compute_lrmsd_fast(method='svd', name=atoms)
        self.lrmsd = sim.compute_lrmsd_pdb2sql(exportpath=None, method='svd', name = atoms)

        # remove intermediate files
        os.system('rm %s/decoy.pdb %s/ref.pdb' %(self.output_dir, self.output_dir))

    def calc_Core_LRMSD(self, reference_pdb, atoms = ['C', 'CA', 'N', 'O']):
        ''' Calculate the L-RMSD between the decoy and reference structure (ground truth)

        Args:
            template_pdb: Bio.PDB object

        Returns: (float) L-RMSD
        '''

        # load target pdb
        if isinstance(reference_pdb, str):  # if its a string, it should be the path of the pdb, then load pdb first
            ref = PDBParser(QUIET=True).get_structure('MHC', reference_pdb)
        else:
            ref = template_pdb

        # pdb2sql needs 1 big chain and 1 ligand chain with correct numbering, for MHCII, this means merging the chains.
        homogenize_pdbs(self.pdb, ref, self.output_dir, anchors = self.target.anchors)

        # Calculate l-rmsd between decoy and reference with pdb2sql
        sim = StructureSimilarity('%s/decoy.pdb' % (self.output_dir.replace(' ', '\\ ')), '%s/ref.pdb' % (self.output_dir))
        self.core_lrmsd = sim.compute_lrmsd_fast(method='svd', name=atoms)

        # remove intermediate files
        os.system('rm %s/decoy.pdb %s/ref.pdb' %(self.output_dir, self.output_dir))


def merge_chains(pdb):
    ''' Merges two chains of MHCII to one chain. pdb2sql can only calculate L-rmsd with one chain.

    Args:
        pdb: Bio.PDB object

    Returns: Bio.PDB object with its M and N chain merged as M chain

    '''
    # Merge chains
    if 'N' in [chain.id for chain in pdb.get_chains()]:

        for j in pdb[0]['N'].get_residues():
            j.id = (j.id[0], j.id[1], 'M')
            pdb[0]['M'].add(j)

        for i in pdb.get_chains():
            for model in pdb:
                for chain in model:
                    if chain.id in ['N']:
                        model.detach_child(chain.id)
    return pdb


def renumber(pdb):
    ''' Renumbers the pdb. Each chain starts at 1

    Args:
        pdb: Bio.PDb object

    Returns:Bio.PDb object with renumbered residues

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


def homogenize_pdbs(decoy, ref, output_dir, anchors =False ):
    ''' Make sure that the decoy and reference structure have the same structure sequences.

    Args:
        decoy: Bio.PDB object of the decoy structure
        ref: Bio.PDB object of the reference structure
        output_dir: (string) directory that is used to write intermediate files

    Returns: (tuple) Bio.PDB objects with the same structure sequence

    '''

    # If you give the anchors, the core L-RMSD will be calculated.
    # The peptide residues before and after the first and last anchor residue will be discarded.
    if anchors:
        for x in range(len(decoy[0]['P'])):
            for i in decoy[0]['P']:
                if i.id[1] < anchors[0] or i.id[1] > anchors[-1]:
                    decoy[0]['P'].detach_child(i.id)
            for i in ref[0]['P']:
                if i.id[1] < anchors[0] or i.id[1] > anchors[-1]:
                    ref[0]['P'].detach_child(i.id)

    # merge chains of the decoy
    decoy = merge_chains(decoy)
    decoy = renumber(decoy)
    # merge chains of the reference
    ref = merge_chains(ref)
    ref = renumber(ref)

    # Write pdbs
    io = PDBIO()
    io.set_structure(decoy)
    io.save('%s/decoy.pdb' % (output_dir))

    io = PDBIO()
    io.set_structure(ref)
    io.save('%s/ref.pdb' % (output_dir))

    return decoy, ref


#
# decoy_path = '/Users/derek/Dropbox/Master_Bioinformatics/Internship/PANDORA_remaster/PANDORA/PANDORA_files/data/outputs/1DLH_1FYT/1FYT.BL00010001.pdb'
# ref_path = '/Users/derek/Dropbox/Master_Bioinformatics/Internship/PANDORA_remaster/PANDORA/PANDORA_files/data/PDBs/pMHCII/1FYT.pdb'
#
# decoy = PDBParser(QUIET=True).get_structure('Decoy', decoy_path)
# ref = PDBParser(QUIET=True).get_structure('Ref', ref_path)
#
#
# decoy_path = '/Users/derek/Dropbox/Master_Bioinformatics/Internship/PANDORA_remaster/PANDORA/PANDORA_files/data/outputs/5KSU_6U3O/6U3O.BL00010001.pdb'
# ref_path = '/Users/derek/Dropbox/Master_Bioinformatics/Internship/PANDORA_remaster/PANDORA/PANDORA_files/data/PDBs/pMHCII/6U3O.pdb'
#
# m = Model(mod.target, mod.output_dir, model_path=decoy_path)
# m.calc_LRMSD(ref_path)
# m.calc_Core_LRMSD(ref_path)
# print(m.lrmsd)
# print(m.core_lrmsd)



