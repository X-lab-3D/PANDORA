from modelling_scripts import data_prep

IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/csv_pkl_files/mhcI_structures_IEDB.csv', 42, 43, ';', empty_rows=[0,1])
