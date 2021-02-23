from Bio import SeqIO
import PANDORA
import os

# target = mod.target
# template = mod.template
# output_dir = '/Users/derek/Dropbox/Master Bioinformatics/Internship/PANDORA_remaster/PANDORA/PANDORA_files/data/outputs/1D9K_1IAK'


class Align:

    def __init__(self, target, template, output_dir = PANDORA.PANDORA_data + '/outputs'):
        '''Performs a alignment of the target and template(s). Will spit out a filename that will be used for modeller.

        :param target: (Target object) The target object that will be aligned to the template structures
        :param template: (Template object) can be either a single Target object or a list of Target objects
        :param output_dir: (string)
        '''
        # Assign target and template. If the template is not given as a list, put in in a list.
        self.target = target
        if isinstance(template, list):
            self.template = template
        else:
            self.template = [template]

        # Create an output file if there isn't one yet: template.id_target.id (chains together the name of the multiple
        # templates if they are given)
        self.output_dir = output_dir + '/%s_%s' %('_'.join([i.PDB_id for i in self.template]), self.target.PDB_id)
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # Find out if the target is MHC class I or II
        self.__MHC_class = target.MHC_class

        # Store the target and template ids for later use
        self.__tar_id = target.PDB_id
        self.__tem_id = [i.PDB_id for i in self.template]

        # Define template m, n and p seqs
        self.__tem_m = [i.chain_seq[0] for i in self.template]
        if self.__MHC_class == 'II':
            self.__tem_n = [i.chain_seq[1] for i in self.template]
        self.__tem_p = [i.peptide for i in self.template]

        # Define target m, n and p seqs. If there are no m and n chains supplied, just take the seqs from the template
        # Only if the target is MHC class II, the n chains are defined.
        if self.target.chain_seq == []:
            self.__tar_m = self.__tem_m[0]
            if self.__MHC_class == 'II':
                self.__tar_n = self.__tem_n[0]
        else:
            self.__tar_m = self.target.chain_seq[0]
            if self.__MHC_class == 'II':
                self.__tar_n = self.target.chain_seq[1]
        self.__tar_p = self.target.peptide



        self.aligned_seqs_and_pept = self.__align_chains()

        self.alignment_file = self.__write_ali_file()


    def __align_chains(self):
        ''' Alignes the chains of MHCI/MHCII structues given a target and (multiple) templates.

        :return: (dict) all aligned chains per structure id {id M: NNNNNNNNNNNN, id N: NNNNNNNNNNNN, id P: NNNNNN}
        '''

        # Align M and N chain for MHC II. Because the target chains need to be aligned to the respective chain of
        # the template, M and N are done seperately and later added together
        if self.__MHC_class == 'II':
            # Align the M chain
            # First write a fasta file containing all chains
            with open('%s/%s_M.fasta' % (self.output_dir,self.__tar_id),"w") as f:
                for i in range(len(self.__tem_id)):
                    f.write('>'+self.__tem_id[i] + ' M\n' + self.__tem_m[i] +'\n')
                f.write('>' + self.__tar_id + ' M\n' + self.__tar_m)
            # Perform MSA with muscle
            msl = MuscleCommandline(input='%s/%s_M.fasta' % (self.output_dir,self.__tar_id), out='%s/%s_M.afa' % (self.output_dir, self.__tar_id))
            os.system(str(msl) + ' -quiet')

            # Align the N chain
            # First write a fasta file containing all chains
            with open('%s/%s_N.fasta' % (self.output_dir,self.__tar_id),"w") as f:
                for i in range(len(self.__tem_id)):
                    f.write('>'+self.__tem_id[i] + ' N\n' + self.__tem_n[i] +'\n')
                f.write('>' + self.__tar_id + ' N\n' + self.__tar_n)
            # Perform MSA with muscle
            msl = MuscleCommandline(input='%s/%s_N.fasta' % (self.output_dir, self.__tar_id), out='%s/%s_N.afa' % (self.output_dir, self.__tar_id))
            os.system(str(msl) + ' -quiet')

            # Merge M and N chain into one file
            os.system('cat %s/*.afa > %s/alignment.afa' %(self.output_dir.replace(' ', '\\ '), self.output_dir.replace(' ', '\\ ')))

        if self.__MHC_class == 'I':
            # Align the M chain
            # First write a fasta file containing all chains
            with open('%s/%s_M.fasta' % (self.output_dir,self.__tar_id),"w") as f:
                for i in range(len(self.__tem_id)):
                    f.write('>'+self.__tem_id[i] + ' M\n' + self.__tem_m[i] +'\n')
                f.write('>' + self.__tar_id + ' M\n' + self.__tar_m)
            # Perform MSA with muscle
            msl = MuscleCommandline(input='%s/%s_M.fasta' % (self.output_dir,self.__tar_id), out='%s/alignment.afa' % (self.output_dir))
            os.system(str(msl) + ' -quiet')

        # Load aligned seqs into dict
        seqs = {v.description: str(v.seq) for (v) in SeqIO.parse('%s/alignment.afa' % (self.output_dir), "fasta")}

        # Align peptides
        # aligned_pepts = {'1D9K P': 'GNSHRGAIEWEGIESG', '1IAK P': '-STDYGILQINSRW--'}
        aligned_pepts = self.align_peptides()

        # Remove all intermediate files
        os.system('rm %s/*.fasta %s/*.afa' % (self.output_dir.replace(' ', '\\ '), self.output_dir.replace(' ', '\\ ')))

        return {**seqs, **aligned_pepts}

    def align_peptides(self):
        ''' align two MHCII peptides based on the 2nd and 3rd anchor position.

        :return: (dict) {'id P': 'aligned_peptide_seq'}
        '''

        if self.__MHC_class == 'II':

            # Make a dict containing {id, (peptide_sequence, [anchors])}
            id_pept_anch = {self.__tar_id: (self.__tar_p, self.target.anchors)}
            for i in self.template:
                id_pept_anch[i.PDB_id] = (i.peptide, i.anchors)

            # Find the peptide with the longest 2nd anchor.
            longest_2nd_anch = max((v[1][1], k) for k, v in id_pept_anch.items())[0]
            # add padding to the left of all shorter peptides, so the 2nd anchor of all peptide aligns
            for k, v in id_pept_anch.items():
                padding_len = longest_2nd_anch - v[1][1]
                id_pept_anch[k] = ('-' * padding_len + v[0], [i + padding_len for i in v[1]])

            # Find the longest 3rd anchor
            # add padding between the 2nd and 3rd anchor if needed to align these
            longest_3nd_anch = max((v[1][2], k) for k, v in id_pept_anch.items())[0]
            for k, v in id_pept_anch.items():
                padding_len = longest_3nd_anch - v[1][2]
                id_pept_anch[k] = (v[0][:v[1][2] - 1] + '-' * padding_len + v[0][v[1][2] - 1:],
                                   v[1][:2] + [i + padding_len for i in v[1][2:]])

            # Add padding to the right of the shortest peptides to match their length
            longest_peptide = max((len(v[0]), k) for k, v in id_pept_anch.items())[0]
            for k, v in id_pept_anch.items():
                padding_len = longest_peptide - len(v[0])
                id_pept_anch[k] = (v[0] + '-' * padding_len, v[1])

            # Reformat
            return {k + ' P': v[0] for k, v in id_pept_anch.items()}


    def __write_ali_file(self):
        ''' Writes the alignement file from self.aligned_seqs_and_pept '''

        seqs = self.aligned_seqs_and_pept

        aligned_seqs = {} #dict for keeping all info that will be written to the .ali file
        for id in set([i.rsplit(' ', 1)[0] for i in seqs.keys()]): #Go through all structures

            # The header
            head = '>P1;'+id

            # The comment line that contains the path to the file, start chain and peptide length for modeller.
            if id in self.__tem_id:
                comment = 'structure:%s:1:M:%s:P::::' % (
                    template[self.__tem_id.index(id)].pdb_path.replace(' ', '\\ '), len(seqs[id + ' P']))
            else:
                comment = 'sequence:::::::::'

            # The sequences
            if self.__MHC_class == 'II':
                seq = '%s/%s/%s*' %(seqs[id+' M'], seqs[id+' N'], seqs[id+' P'])
            if self.__MHC_class == 'I':
                seq = '%s/%s*' % (seqs[id + ' M'], seqs[id + ' P'])

            aligned_seqs[id] = (head, comment, seq)

        # Write actual .ali file
        alignment_file = '%s/%s.ali' %(self.output_dir, self.target.PDB_id)
        with open(alignment_file, 'w') as f:
            for k,v in aligned_seqs.items():
                f.write(v[0]+'\n'+v[1]+'\n'+v[2]+'\n\n')

        return alignment_file


# A = Align(mod.target, mod.template)








#
#     return ali_p_target, ali_p_template