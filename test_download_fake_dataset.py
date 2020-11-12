from modelling_scripts import data_prep

IDD, bad_IDs = data_prep.imgt_retrieve_clean('data/csv_pkl_files/fake_dataset_5TEZ.csv', 42, 43, ';', empty_rows=[0,1])
