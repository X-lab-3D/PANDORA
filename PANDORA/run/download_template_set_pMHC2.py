import sys
sys.path.append('./')

from PANDORA.parsing import url_protocols
from PANDORA.parsing import structures_parser


#data_prep.download_ids_alleles_imgt(out_tsv = 'auto_generated_IDs_alleles_from_IMGT.tsv', out_pkl = 'IDs_and_alleles_identity_percs_from_imgt.pkl', print_outfiles = False)
#IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/csv_pkl_files/auto_generated_IDs_alleles_from_IMGT.tsv', 0, 1, '\t')

print()
print('###################################')
print('Downloading structure IDs from IMGT')
print('###################################')
print()

IDs_list = url_protocols.download_ids_imgt('MH2', out_tsv='all_MH2_IDs.tsv')


print()
print('###################################')
print('Parsing PDB files')
print('###################################')
print()

IDs_dict, bad_IDs = structures_parser.parse_pMHCII_pdbs(IDs_list)
