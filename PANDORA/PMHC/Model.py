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
            self.pdb = PDBParser(QUIET=True).get_structure(self.target.PDB_id, self.model_path)


    def calc_LRMSD(self, reference_pdb):
        ''' Calculate the L-RMSD between the decoy and reference structure (ground truth)

        Args:
            template_pdb: Bio.PDB object

        Returns: (float) L-RMSD
        '''

        # load target pdb
        if isinstance(reference_pdb, str):  # if its a string, it should be the path of the pdb, then load pdb first
            ref = PDBParser(QUIET=True).get_structure('MHC', reference_pdb)
        else:            ref = template_pdb

        # make sure the structure sequences of the decoy and template are the same. If not, make them the same by
        # removing some residues. This can happen if there is no sequence if the target structure given by the user and
        # the template structure is used.
        homogenize_pdbs(self.pdb, ref, self.output_dir)

        # Calculate l-rmsd between decoy and reference with pdb2sql
        sim = StructureSimilarity('%s/decoy.pdb' % (self.output_dir.replace(' ', '\\ ')), '%s/ref.pdb' % (self.output_dir))
        self.lrmsd = sim.compute_lrmsd_pdb2sql(exportpath=None,
                                          method='svd')
        # remove intermediate files
        os.system('rm %s/decoy.pdb %s/ref.pdb' %(self.output_dir, self.output_dir))


def merge_chains(pdb):
    ''' Merges two chains of MHCII to one chain. pdb2sql can only calculate L-rmsd with one chain.

    Args:
        pdb: Bio.PDB object

    Returns: Bio.PDB object with its M and N chain merged as M chain

    '''
    # Merge chains
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


def homogenize_pdbs(decoy, ref, output_dir):
    ''' Make sure that the decoy and reference structure have the same structure sequences.

    Args:
        decoy: Bio.PDB object of the decoy structure
        ref: Bio.PDB object of the reference structure
        output_dir: (string) directory that is used to write intermediate files

    Returns: (tuple) Bio.PDB objects with the same structure sequence

    '''


    with open('%s/decoy_tar.fasta' % (output_dir), "w") as f:
        f.write('>decoy\n' + seq1(''.join([i.resname for i in decoy[0]['M'].get_residues()])) + seq1(''.join([i.resname for i in decoy[0]['N'].get_residues()])) + '\n')
        f.write('>ref\n' + seq1(''.join([i.resname for i in ref[0]['M'].get_residues()])) + seq1(''.join([i.resname for i in ref[0]['N'].get_residues()])) + '\n')

    # Perform MSA with muscle
    msl = MuscleCommandline(input='%s/decoy_tar.fasta' % (output_dir),
                            out='%s/decoy_tar.afa' % (output_dir))
    os.system(str(msl) + ' -quiet')

    seqs = {v.description: str(v.seq) for (v) in SeqIO.parse('%s/decoy_tar.afa' % (output_dir), "fasta")}

    os.system('rm %s/decoy_tar.afa %s/decoy_tar.fasta' % (output_dir,output_dir))

    # merge chains of the decoy
    decoy = merge_chains(decoy)
    # Remove residues based on the alignment
    res_to_remove = [i for i in range(0, len(seqs['ref'])) if seqs['ref'][i] == '-']
    resnr = 0
    for res in decoy[0]['M']:
        # if seq1(res.resname) != seqs['decoy'][resnr]:
        if resnr in res_to_remove:
            decoy[0]['M'].detach_child(res.id)
        resnr += 1
    # Renumber
    decoy = renumber(decoy)

    # merge chains of the reference
    ref = merge_chains(ref)
    # Remove residues based on the alignment
    res_to_remove = [i for i in range(0, len(seqs['decoy'])) if seqs['decoy'][i] == '-']
    resnr = 0
    for res in ref[0]['M']:
        # if seq1(res.resname) != seqs['decoy'][resnr]:
        if resnr in res_to_remove:
            ref[0]['M'].detach_child(res.id)
        resnr += 1
    # renumber
    ref = renumber(ref)

    # Write pdbs
    io = PDBIO()
    io.set_structure(decoy)
    io.save('%s/decoy.pdb' % (output_dir))

    io = PDBIO()
    io.set_structure(ref)
    io.save('%s/ref.pdb' % (output_dir))

    return decoy, ref



# for i in range(len(logf)):
#
#     decoy = PDBParser(QUIET=True).get_structure('MHC', mod.output_dir + '/' + logf[i][0])
#     ref = PDBParser(QUIET=True).get_structure('MHC', '/Users/derek/Dropbox/Master Bioinformatics/Internship/PANDORA_remaster/PANDORA/PANDORA_files/data/PDBs/pMHCII/1IAK.pdb')
#
#     homogenize_pdbs(decoy, ref, mod.output_dir)
#
#     sim = StructureSimilarity('%s/decoy.pdb' % (mod.output_dir),'%s/ref.pdb' % (mod.output_dir))
#
#     lrmsd = sim.compute_lrmsd_pdb2sql(exportpath=None,
#         method='svd')
#
#     print(lrmsd)
#
#
#
#
# logf = []
# f = open(mod.output_dir + '/modeller.log')
# for line in f:
#     if line.startswith('1IAK.'):
#         l = line.split()
#         if len(l) > 2:
#             logf.append(tuple(l))
# f.close()
#
# henk = Model(mod.target, mod.output_dir + '/' + logf[0][0], molpdf=logf[0][1], dope=logf[0][2])
#
# henk.pdb

# for i in logf:
#     decoy = mod.output_dir + '/' + i[0]
#     ref = '/Users/derek/Dropbox/Master Bioinformatics/Internship/PANDORA_remaster/PANDORA/PANDORA_files/data/PDBs/pMHCII/1IAK.pdb'
#
#     homogenize_pdbs(decoy, ref, mod.output_dir)
#
#     sim = StructureSimilarity('%s/decoy.pdb' % (mod.output_dir), '%s/ref.pdb' % (mod.output_dir))
#
#     lrmsd = sim.compute_lrmsd_pdb2sql(exportpath=None,
#                                       method='svd')
#     print(lrmsd)
#
#
#
