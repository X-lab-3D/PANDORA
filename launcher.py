###    ###
from scripts import data_prep
import pickle
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
import time
from random import choice

#IDD = data_prep.imgt_retrieve_clean('data/final_mhc1_3d_structure_data_with_pdb_ids.tsv')

## TODO:
## Retrieve IDs : Alleles
## Read user's input (seq and allele)
## Select template structure and get its sequence
## Compare User's input with same allele sequences (select with identity the best one)


IDD_file = open('data/IDs_ChainsCounts_dict.pkl', 'rb')
IDD = pickle.load(IDD_file)
IDD_file.close()

allele_ID = {}
for key in IDD:
    try:
        allele_ID[IDD[key]['allele']].append(key)
    except KeyError:
        allele_ID[IDD[key]['allele']] = [key]

print('##############################')
print('Please enter your allele name. You can choose between the ones listed below')
print('##############################')
print('')
print('')
print(allele_ID.keys())
#allele = input()
#allele, inseq_file = 'FLA-E*01801', 'data/5xmf.fasta'
allele, inseq_file = 'H2-Kb', 'data/3roo.fasta'

print('##############################')
print('Please select your sequence file.')
print('##############################')
print('')

#inseq_file = 'data/5xmf.fasta'
#inseq_file = 'data/1tom.fasta'
inseqs = list(SeqIO.parse(inseq_file, 'fasta'))
A_target = inseqs[0]
P_target = inseqs[1]

sequences = data_prep.get_pdb_seq(allele_ID[allele])
to_write = []
for i, ID in enumerate(allele_ID[allele]):
    template = SeqRecord(Seq(sequences[i][0], IUPAC.protein), id=ID, name = allele + ID)
    to_write.append(template)
    SeqIO.write((template, A_target), "data/FASTAs/%s.fasta" %ID, "fasta")


if len(allele_ID[allele]) == 1:
    template_ID = allele_ID[allele][0]
    #subprocess.check_call(['muscle', '-in %s.fasta -clw' %template_ID])
    t1 = time.time()
    os.system('muscle -in data/FASTAs/%s.fasta -out data/Alignments/try.afa' %ID)
    tf = time.time() - t1
    print('IT TOOK: ')
    print(tf)
else:
    score_dict = {}
    max_id = 0
    putative_templates = []
    for ID in allele_ID[allele]:
        alifile = 'data/Alignments/%s.afa' %ID
        muscle_commands = ['muscle', '-in', 'data/FASTAs/%s.fasta' %ID, '-out', '%s' %alifile]
        #os.system('muscle -in data/FASTAs/%s.fasta -out %s' %(ID, alifile))
        subprocess.check_call(muscle_commands)
        os.system('rm data/FASTAs/%s.fasta' %ID)
        aligned = list(SeqIO.parse(alifile, 'fasta'))
        id_count = 0
        for x, y in zip((aligned[0].seq), str(aligned[1].seq)):
            if x == y:
                id_count += 1
            else:
                pass
        if id_count > max_id:
            max_id = id_count
        score_dict[ID] = id_count
        print('ID COUNT:  ', id_count)
    for key in score_dict:
        if score_dict[key] == max_id:
            putative_templates.append(key)
        else:
            os.system('rm data/Alignments/%s.afa' %key)
    template_ID = choice(list(score_dict))
            
        