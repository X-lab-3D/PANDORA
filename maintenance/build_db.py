## import requested modules
from PANDORA import Database


## Create local Database

db_path = '~/PANDORA_databases/default'
Database.create_db_folders(db_path=db_path)

db = Database.Database()
db.construct_database(save=f'{db_path}/database/PANDORA_database.pkl', 
                        data_dir = db_path,            
                        n_jobs=64, download=True)
