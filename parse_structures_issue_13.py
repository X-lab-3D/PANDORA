from modelling_scripts import url_protocols
from modelling_scripts import structures_parser
import csv

#data_prep.download_ids_alleles_imgt(out_tsv = 'auto_generated_IDs_alleles_from_IMGT.tsv', out_pkl = 'IDs_and_alleles_identity_percs_from_imgt.pkl', print_outfiles = False)
#IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/csv_pkl_files/auto_generated_IDs_alleles_from_IMGT.tsv', 0, 1, '\t')

ids_filename = 'data/csv_pkl_files/mhcI_structures_IEDB.csv'
id_clmn = 42
empty_rows = [0,1]
delimiter = ';'

IDs = []
with open(ids_filename, 'r') as idsfile:
    spamreader = csv.reader(idsfile, delimiter=delimiter)
    for i, row in enumerate(spamreader):
        if i in empty_rows:
            pass
        else:
            ID = row[id_clmn]
            IDs.append(ID)


print()
print('###################################')
print('Parsing PDB files')
print('###################################')
print()

IDs_dict, bad_IDs = structures_parser.parse_pMHCI_pdbs(IDs)
