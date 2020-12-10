from modelling_scripts import data_prep

#data_prep.download_ids_alleles_imgt(out_tsv = 'auto_generated_IDs_alleles_from_IMGT.tsv', out_pkl = 'IDs_and_alleles_identity_percs_from_imgt.pkl', print_outfiles = False)
IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/csv_pkl_files/auto_generated_IDs_alleles_from_IMGT.tsv', 0, 1, '\t')
