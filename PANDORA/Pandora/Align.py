from Bio import SeqIO
from Bio.Align import substitution_matrices
from Bio.Align import PairwiseAligner
import Bio.Align
import PANDORA
import os
import subprocess


class Align:

    def __init__(self, target, template, clip_C_domain=False,
                 output_dir=PANDORA.PANDORA_data + '/outputs', remove_terms=True):
        ''' Performs a alignment of the target and template(s). Returns a filename that will be used for modeller.

        Args:
            target: (Target object) The target object that will be aligned to
                the template structures
            template: (list): list with a single Template object
            clip (bool or list): if True, clips away the C-like domain, levaing only
                the G-domain according to IMGT. If a listcontaining the G domain(s) 
                span is provided, will use it to cut the sequence. The list should have 
                this format: [(1,182)] for MHCI and [(1,91),(1,86)] for MHCII.
            output_dir: (string) output directory
            remove_terms (bool): if True, removes N-terminal and C-terminal
                regions not present in the template sequence
        '''

        self.clip_C_domain = clip_C_domain
        # Assign target and template. If the template is not given as a list, put in in a list.
        self.target = target
        if type(template) == list:
            self.template = template
        else:
            self.template = [template]
        self.output_dir = output_dir # + '/%s_%s' %('_'.join([i.id for i in self.template]), self.target.id)
        self.__muscle_command__ = 'muscle -align %s -output %s -quiet'


        # Find out if the target is MHC class I or II
        self.MHC_class = target.MHC_class

        # Store the target and template ids for later use
        self.tar_id = target.id
        self.tem_id = [i.id for i in self.template]

        # Define template m, n and p seqs
        self.tem_m = self.template[0].M_chain_seq
        if self.MHC_class == 'II':
            self.tem_n = self.template[0].N_chain_seq
        self.tem_p = self.template[0].peptide 
        
        # if clip, cut away the C-like domain(s)
        if self.clip_C_domain == False:
            pass
        
        elif self.clip_C_domain == True:
            #try:
            self.tem_m, self.tem_m_c = self.clip(self.tem_m, self.template[0].G_domain_span[0])
            if self.MHC_class == 'II':
                self.tem_n, self.tem_n_c = self.clip(self.tem_n, self.template[0].G_domain_span[1])
            #except:
            #    print('WARNING: an error occurred while clipping the template sequence.')
                
        elif type(self.clip_C_domain) == list:
            #try:
            self.tem_m, self.tem_m_c = self.clip(self.tem_m, self.clip_C_domain[0])
            if self.MHC_class == 'II':
                self.tem_n, self.tem_n_c = self.clip(self.tem_n, self.clip_C_domain[1])
            #except:
            #    print('WARNING: an error occurred while clipping the template sequence.')

        # Define target m, n and p seqs. If there are no m and n chains supplied, just take the seqs from the template
        if self.target.M_chain_seq == '':
            print('WARNING: No M chain sequence could be retrieved for target %s' %self.target.id)
            print('PANDORA will use the M chain sequence from the best template.')
            print('To avoid this, please provide an M chain sequence to your target.')
            self.tar_m = self.tem_m
        else:
            self.tar_m = self.target.M_chain_seq

        # N chains (only for MHCII)
        if self.MHC_class == 'II' and self.target.N_chain_seq == '':
            print('WARNING: No M chain sequence could be retrieved for target %s' %self.target.id)
            print('PANDORA will use the M chain sequence from the best template.')
            print('To avoid this, please provide an M chain sequence to your target.')
            self.tar_n = self.tem_n
        elif self.MHC_class == 'II' and self.target.N_chain_seq != '':
            self.tar_n = self.target.N_chain_seq

        # P chain
        self.tar_p = self.target.peptide

        # Perform alignment
        self.aligned_seqs_and_pept = self.align_chains()
        # Cut extra N-terminal and C-terminal
        if self.clip_C_domain==True or remove_terms==True:
            self.remove_terms()
        if self.clip_C_domain:
            self.aligned_seqs_and_pept['%s M' %self.template[0].id] += self.tem_m_c
            self.aligned_seqs_and_pept['%s M' %self.target.id] += ('').join(
                ['-' for i in range(len(self.tem_m_c))])
            if self.MHC_class == 'II':
                self.aligned_seqs_and_pept['%s N' %self.template[0].id] += self.tem_n_c
                self.aligned_seqs_and_pept['%s N' %self.target.id] += ('').join(
                ['-' for i in range(len(self.tem_n_c))])
        #Write alignment file for MODELLER
        self.alignment_file = self.write_ali_file()

    def clip(self, seq, span):
        seq_1 = seq[span[0]-1:span[1]]
        seq_2 = seq[span[1]:]
        return seq_1, seq_2

    def align_chains(self):
        ''' Alignes the chains of MHCI/MHCII structues given a target and (multiple) templates.

        Returns: (dict) all aligned chains per structure id {id M: NNNNNNNNNNNN, id N: NNNNNNNNNNNN, id P: NNNNNN}

        '''

        # Align M and N chain for MHC II. Because the target chains need to be aligned to the respective chain of
        # the template, M and N are done seperately and later added together
        if self.MHC_class == 'I':
            chains = {"M" : 'alignment'}
        elif self.MHC_class == 'II':
            chains = {"M":f'{self.tar_id}_M', 
                    "N":f'{self.tar_id}_N'}

        for chain, afa_name in chains.items():
            # Align the M chain
            # First write a fasta file containing all chains
            with open(f'{self.output_dir}/{self.tar_id}_{chain}.fasta',"w") as f:
                for i in range(len(self.tem_id)):
                    # Write template id \n template seq
                    f.write(f'>{self.tem_id[i]} {chain}\n{self.tem_m}\n')
                # Write target id \n target seq
                f.write(f'>{self.tar_id} {chain}\n{self.tar_m}')
            # Perform MSA with muscle
            in_file_muscle = f'{self.output_dir}/{self.tar_id}_{chain}.fasta'
            out_file_muscle = f'{self.output_dir}/{afa_name}.afa'
            p = subprocess.check_call(self.__muscle_command__ % (in_file_muscle,out_file_muscle),shell=True)

        if self.MHC_class == 'II':
            # Merge M and N chain into one file
            os.system('cat %s/%s_M.afa %s/%s_N.afa > %s/alignment.afa' % (
            self.output_dir.replace(' ', '\\ '), self.tar_id, self.output_dir.replace(' ', '\\ '), self.tar_id,
            self.output_dir.replace(' ', '\\ ')))

        # Load aligned seqs into dict
        seqs = {v.description: str(v.seq) for (v) in SeqIO.parse('%s/alignment.afa' % (self.output_dir), "fasta")}

        # Align peptides
        # aligned_pepts = {'1D9K P': 'GNSHRGAIEWEGIESG', '1IAK P': '-STDYGILQINSRW--'}
        aligned_pepts = self.align_peptides()

        # Remove all intermediate files
        os.system('rm %s/*.fasta %s/*.afa' % (self.output_dir.replace(' ', '\\ '), self.output_dir.replace(' ', '\\ ')))
        return {**seqs, **aligned_pepts}

    def align_peptides(self):
        ''' Align MHCI peptides based on their anchors. This function adds padding to the left and the right of
            the anchors and uses muscle to align the core. For MHCII, the peptides are aligned on the 2nd and 3rd anchor
            position. These two anchor positions are present in all MHCII structures (1 and 4 are missing sometimes).
            For MHCII, no Muscle alignment is performed because there is very little distance between the 2nd and 3rd
            anchor.

        Returns: (dict) {'id P': 'aligned_peptide_seq'}

        '''

        # Make a dict containing {id, (peptide_sequence, [anchors])}
        id_pept_anch = {self.tar_id: (self.tar_p, self.target.anchors)}
        for i in self.template:
            id_pept_anch[i.id] = (i.peptide, i.anchors)

        if self.MHC_class == 'I':

            # Find the peptide with the longest 1st anchor.
            longest_1st_anch = max((v[1][0], k) for k, v in id_pept_anch.items())[0]
            # add padding to the left of all shorter peptides, so the 1st anchor of all peptide aligns
            for k, v in id_pept_anch.items():
                padding_len = longest_1st_anch - v[1][0]
                id_pept_anch[k] = ('-' * padding_len + v[0], [i + padding_len for i in v[1]])

            # align peptide cores
            cores = {k: v[0][v[1][0]:(v[1][1]) - 1] for k, v in id_pept_anch.items()}

            # Check if all cores are of equal length, otherwise align them
            if len(set([len(v) for k,v in cores.items()])) != 1:

                # Write Fasta file
                with open('%s/pept_cores.fasta' % (self.output_dir), "w") as f:
                    for k, v in cores.items():
                        f.write('>' + k + '\n' + v + '\n')

                # Run Muscle
                in_file_muscle = '%s/pept_cores.fasta' % (self.output_dir)
                out_file_muscle = '%s/pept_cores.afa' % (self.output_dir)
                p = subprocess.call(self.__muscle_command__ % (in_file_muscle,out_file_muscle),shell=True)

                # Load aligned seqs into dict
                aligned_cores = {v.description: str(v.seq) for (v) in
                                 SeqIO.parse('%s/pept_cores.afa' % (self.output_dir), "fasta")}

                # Remove intermediate files
                # os.system(
                #     'rm %s/*.fasta %s/*.afa' % (
                #     self.output_dir.replace(' ', '\\ '), self.output_dir.replace(' ', '\\ ')))

                # Add aligned cores in between anchors
                for k, v in id_pept_anch.items():
                    id_pept_anch[k] = (v[0][:v[1][0]] + aligned_cores[k] + v[0][(v[1][1]) - 1:], v[1])


        if self.MHC_class == 'II':

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


    def remove_terms(self):
        '''Removes N-terminal and C-terminal regions not present in the template'''
        if self.MHC_class == 'I':
            chains = ['M']
        elif self.MHC_class == 'II':
            chains = ['M', 'N']

        
        def cut_excess(ref_seq, ref_id, decoy_id, chain):
            
            #N terminal excess
            N_term_excess = 0
            if ref_seq.startswith('-'):
                for aa in ref_seq:
                    if aa == '-':
                        N_term_excess += 1
                    else:
                        break
    
                #Cut N-terminal excess regions
                self.aligned_seqs_and_pept['%s %s' %(ref_id, chain)] = self.aligned_seqs_and_pept['%s %s' %(ref_id, chain)][N_term_excess:]
                self.aligned_seqs_and_pept['%s %s' %(decoy_id, chain)] = self.aligned_seqs_and_pept['%s %s' %(decoy_id, chain)][N_term_excess:]
    
    
            #C terminal excess
            C_term_excess = 0
            if ref_seq.endswith('-'):
                for aa in templ_seq[::-1]:
                    if aa == '-':
                        C_term_excess += 1
                    else:
                        break
    
                #Cut C-terminal excess regions
                self.aligned_seqs_and_pept['%s %s' %(ref_id, chain)] = self.aligned_seqs_and_pept['%s %s' %(ref_id, chain)][:-C_term_excess]
                self.aligned_seqs_and_pept['%s %s' %(decoy_id, chain)] = self.aligned_seqs_and_pept['%s %s' %(decoy_id, chain)][:-C_term_excess]
            print('Chain %s, excess N %i, C %i' %(chain, N_term_excess, C_term_excess))
        
        for chain in chains:
            templ_id = self.template[0].id
            target_id = self.target.id
            templ_seq = self.aligned_seqs_and_pept['%s %s' %(templ_id, chain)]
            #target_seq = self.aligned_seqs_and_pept['%s %s' %(target_id, chain)]
            cut_excess(ref_seq=templ_seq, ref_id=templ_id, decoy_id=target_id, chain=chain)
            #cut_excess(ref_seq=target_seq, ref_id=target_id, decoy_id=templ_id, chain=chain)

    def write_ali_file(self):
        ''' Writes the alignement file from self.aligned_seqs_and_pept '''

        seqs = self.aligned_seqs_and_pept

        aligned_seqs = {} #dict for keeping all info that will be written to the .ali file
        for id in set([i.rsplit(' ', 1)[0] for i in seqs.keys()]): #Go through all structures

            # The header
            head = '>P1;'+id

            # The comment line that contains the path to the file, start chain and peptide length for modeller.
            if id in self.tem_id:
                comment = 'structure:%s:1:M:%s:P::::' % (
                    os.path.basename(self.template[self.tem_id.index(id)].get_pdb_path()), len(seqs[id + ' P']))
            else:
                comment = 'sequence:::::::::'

            # The sequences
            if self.MHC_class == 'II':
                seq = '%s/%s/%s*' %(seqs[id+' M'], seqs[id+' N'], seqs[id+' P'])
            if self.MHC_class == 'I':
                seq = '%s/%s*' % (seqs[id + ' M'], seqs[id + ' P'])

            aligned_seqs[id] = (head, comment, seq)

        # Write actual .ali file
        alignment_file = '%s/%s.ali' %(self.output_dir, self.target.id)
        with open(alignment_file, 'w') as f:
            for k,v in aligned_seqs.items():
                f.write(v[0]+'\n'+v[1]+'\n'+v[2]+'\n\n')

        return alignment_file


class Align2:

    def init(self, target, template, output_dir=PANDORA.PANDORA_data + '/outputs'):
        ''' Experimental function. Not used. 
            Performs a alignment of the target and template(s). Will spit out a filename that will be used for modeller.

        Args:
            target: (Target object) The target object that will be aligned to the template structures
            template: (Template object) can be either a single Target object or a list of Target objects
            output_dir: (string)
        '''

        self.target = target
        if isinstance(template, list):
            self.templates = template
        else:
            self.templates = [template]
        self.output_dir = output_dir
        self.__muscle_command__ = 'muscle -align %s -output %s -quiet'



    def align_templates(self):
        ''' Alignes the chains of MHCI/MHCII structues given a target and (multiple) templates.

        Returns: (dict) all aligned chains per structure id {id M: NNNNNNNNNNNN, id N: NNNNNNNNNNNN, id P: NNNNNN}

        '''
        # Define some variables
        MHC_class = self.target.MHC_class
        opd = self.output_dir
        tar_id = self.target.id
        tem_id = [i.id for i in self.templates]

        tem_m = [i.M_chain_seq for i in self.templates]
        if MHC_class == 'II':
            tem_n = [i.N_chain_seq for i in self.templates]
        tem_p = [i.peptide for i in self.templates]

        tar_p = self.target.peptide
        if self.target.M_chain_seq == '':
            tar_m = tem_m[0]
            if MHC_class == 'II' and self.target.N_chain_seq == '':
                tar_n = tem_n[0]
        else:
            tar_m = self.target.M_chain_seq
            if MHC_class == 'II':
                tar_n = self.target.N_chain_seq


        # Align M and N chain for MHC II. Because the target chains need to be aligned to the respective chain of
        # the template, M and N are done seperately and later added together
        if self.MHC_class == 'I':
            chains = {"M" : 'alignment'}
        elif self.MHC_class == 'II':
            chains = {"M":'M', "N":'N'}

        for chain, afa_name in chains.items():
            # Align the M chain
            # First write a fasta file containing all chains
            with open(f'{opd}/{self.tar_id}_{chain}.fasta',"w") as f:
                for i in range(len(self.tem_id)):
                    # Write template id \n template seq
                    f.write(f'>{self.tem_id[i]} {chain}\n{self.tem_m}\n')
                # Write target id \n target seq
                f.write(f'>{self.tar_id} {chain}\n{self.tar_m}')
            # Perform MSA with muscle
            in_file_muscle = f'{opd}/{self.tar_id}_{chain}.fasta'
            out_file_muscle = f'{opd}/{self.tar_id}_{afa_name}.afa'
            p = subprocess.check_call(self.__muscle_command__ % (in_file_muscle,out_file_muscle),shell=True)

        if self.MHC_class == 'II':
            # Merge M and N chain into one file
            os.system('cat %s/%s_M.afa %s/%s_N.afa > %s/alignment.afa' % (
            opd.replace(' ', '\\ '), self.tar_id, opd.replace(' ', '\\ '), self.tar_id,
            opd.replace(' ', '\\ ')))

        # Align peptides
        # aligned_pepts = {'1D9K P': 'GNSHRGAIEWEGIESG', '1IAK P': '-STDYGILQINSRW--'}
        self.align_peptides()

        self.aligned_seqs = {**seqs, **self.aligned_pepts}

        # Remove all intermediate files
        os.system('rm %s/*.fasta %s/*.afa' % (opd, opd))

        # Using the aligned sequences, write the alignment file.
        self.alignment_file = self.write_ali_file()


    def align_peptides(self, subt_matrix='PAM30'):
        ''' This function alignes all template peptides with the target peptide based on the anchors (pairwise).
            For MHCI, the left (L), anchors + middle (A,M) and right (R) part of the peptide are aligned and scored:
             LL AMMMMMA RR. This way the anchors always line up. In the end the three parts and their scores are added
             together.
            For MHCII, the same thing is done. However, instead the part between the 2nd and 3rd anchor is taken as as
            middle part. All IMGT structures have proper 2nd and 3rd anchors, while their 1st and especially 4th anchor
            can sometimes be debatable.

        Args:
            subt_matrix: Which subtitution matrix to use? PAM30, BLOSUM80, etc

        Returns: self.aligned_pepts and self.pept_alignment_scores

        '''

        # keep track of alignment scores
        scores = {}
        MHC_class = self.target.MHC_class
        tar_pept = self.target.peptide
        tar_anchors = self.target.anchors
        # Perform a pairwise alignment of the target and all templates for the MHC M chain and peptide
        for i in self.templates + [self.target]:
            temp_anchors = i.anchors

            if MHC_class == 'I':
                # Take the left (of anchor 1) middle (between anchor 1 and 2) and right (of anchor 2) of target peptide
                tar_left = tar_pept[:tar_anchors[0] - 1]
                tar_right = tar_pept[tar_anchors[1]:]
                tar_middle = tar_pept[tar_anchors[0] - 1:tar_anchors[1]]
                # Take the left (of anchor 1) middle (between anchor 1 and 2) and right (of anchor 2) of template peptide
                temp_left = i.peptide[:temp_anchors[0] - 1]
                temp_right = i.peptide[temp_anchors[1]:]
                temp_middle = i.peptide[temp_anchors[0] - 1:temp_anchors[1]]

            if MHC_class == 'II':
                # Take the left (of anchor 2) middle (between anchor 2 and 3) and right (of anchor 3) of target peptide
                tar_left = tar_pept[:tar_anchors[1] - 1]
                tar_right = tar_pept[tar_anchors[2]:]
                tar_middle = tar_pept[tar_anchors[1] - 1:tar_anchors[2]]
                # Take the left (of anchor 2) middle (between anchor 2 and 3) and right (of anchor 3) of template peptide
                temp_left = i.peptide[:temp_anchors[1] - 1]
                temp_right = i.peptide[temp_anchors[2]:]
                temp_middle = i.peptide[temp_anchors[1] - 1:temp_anchors[2]]

            # Create aligner object. The high gap penalties make sure gaps are only created if really needed.
            # For example if the anchors of target and template do not match and a gap needs to be introduced.
            aligner = Bio.Align.PairwiseAligner()
            aligner.substitution_matrix = substitution_matrices.load(subt_matrix)  # PAM30 for pept??
            aligner.gap_score = -1000
            aligner.end_open_gap_score = -1000000
            aligner.internal_open_gap_score = -10000

            # Aligns the middle part of the target and template peptide
            aligned_middle = aligner.align(tar_middle, temp_middle)
            aligned_middle_score = aligned_middle.score
            aligned_middle_seq = format([a for a in aligned_middle][0]).split('\n')[2]

            # Align residues left of the first anchors (if there are any)
            aligned_left_score = 0
            if len(tar_left) > 0 and len(temp_left) > 0:
                # Create aligner object. Make sure no gaps are created in the right side of the sequence
                aligner = Bio.Align.PairwiseAligner()
                aligner.substitution_matrix = substitution_matrices.load(subt_matrix)
                aligner.gap_score = -1000
                aligner.right_open_gap_score = -1000000
                aligner.internal_open_gap_score = -10000
                # Align the sequences
                aligned_left = aligner.align(tar_left, temp_left)
                aligned_left_score = aligned_left.score
                aligned_left_seq = format([a for a in aligned_left][0]).split('\n')[2]
            elif len(temp_left) > 0:
                aligned_left_seq = temp_left
            else:
                aligned_left_seq = ''

            # Align residues right of the second anchors (if there are any)
            aligned_right_score = 0
            if len(tar_right) > 0 and len(temp_right) > 0:
                # Create aligner object. Make sure no gaps are created in the left side of the sequence
                aligner = Bio.Align.PairwiseAligner()
                aligner.substitution_matrix = substitution_matrices.load(subt_matrix)
                aligner.gap_score = -1000
                aligner.left_open_gap_score = -1000000
                aligner.internal_open_gap_score = -10000
                # Align the sequences
                aligned_right = aligner.align(tar_right, temp_right)
                aligned_right_score = aligned_right.score
                aligned_right_seq = format([a for a in aligned_right][0]).split('\n')[2]
            elif len(temp_left) > 0:
                aligned_right_seq = temp_right
            else:
                aligned_right_seq = ''

            # Calculate the total similarity score by summing the scores of the left, middle and right part.
            total_score = aligned_left_score + aligned_middle_score + aligned_right_score
            # Paste the left, middle and right part of the template peptide together again.
            aligned_seq = aligned_left_seq + aligned_middle_seq + aligned_right_seq

            # add results to dict
            scores[i.id] = [total_score, aligned_seq, i.anchors]

        # Add padding left and right to make sure all peptides are the same length
        if MHC_class == 'I':
            scores = self.add_padding(scores, first_aligned_anch_idx=0)
        if MHC_class == "II":
            scores = self.add_padding(scores, first_aligned_anch_idx=1)

        # Dict that is later used in align_templates()
        self.aligned_pepts = {k + ' P': v[1] for k, v in scores.items()}

        # Remove the score of the self alignment between target.peptide and target.peptide. Also remove gaps
        scores.pop(self.target.id)
        # Dict of aligned peptides with their scores
        self.pept_alignment_scores = {k:(v[0],v[1], v[2]) for k,v in scores.items()}


    def add_padding(self, scores, first_aligned_anch_idx=1):
        ''' Internal function for align_peptides() to add --padding-- left and right of peptides. Adds -- according to
            anchors positions

        Args:
            scores: dict: {'1FYT': (-2050.0, '---PKYVKQNTLKLAT--', [3, 6, 8, 11])
            first_aligned_anch_idx: int, 0 for MHCI and 1 for MHCII

        Returns: scores dict with padded peptides

        '''

        # add padding to the left of all shorter peptides, so the 1st (MHCI) or 2nd (MHCII) anchor of all peptide aligns
        longest_2nd_anch = max((v[2][first_aligned_anch_idx], k) for k, v in scores.items())[0]
        for k, v in scores.items():
            padding_len = longest_2nd_anch - v[2][first_aligned_anch_idx]
            scores[k] = (v[0], '-' * padding_len + v[1], [i + padding_len for i in v[2]])

        # Add padding to the right of the shortest peptides to match their length
        longest_peptide = max((len(v[1]), k) for k, v in scores.items())[0]
        for k, v in scores.items():
            padding_len = longest_peptide - len(v[1])
            scores[k] = (v[0], v[1] + '-' * padding_len, v[2])

        return scores

    def write_ali_file(self):
        ''' Writes the alignement file from self.aligned_seqs_and_pept '''

        seqs = self.aligned_seqs
        tem_id = [i.id for i in self.templates]

        aligned_seqs = {} #dict for keeping all info that will be written to the .ali file
        for id in set([i.rsplit(' ', 1)[0] for i in seqs.keys()]): #Go through all structures

            # The header
            head = '>P1;'+id

            # The comment line that contains the path to the file, start chain and peptide length for modeller.
            if id in tem_id:
                comment = 'structure:%s:1:M:%s:P::::' % (
                    os.path.basename(self.templates[tem_id.index(id)].get_pdb_path()), len(seqs[id + ' P']))
            else:
                comment = 'sequence:::::::::'

            # The sequences
            if self.target.MHC_class == 'II':
                seq = '%s/%s/%s*' %(seqs[id+' M'], seqs[id+' N'], seqs[id+' P'])
            if self.target.MHC_class == 'I':
                seq = '%s/%s*' % (seqs[id + ' M'], seqs[id + ' P'])

            aligned_seqs[id] = (head, comment, seq)

        # Write actual .ali file
        alignment_file = '%s/%s.ali' %(self.output_dir, self.target.id)
        with open(alignment_file, 'w') as f:
            for k,v in aligned_seqs.items():
                f.write(v[0]+'\n'+v[1]+'\n'+v[2]+'\n\n')

        return alignment_file
