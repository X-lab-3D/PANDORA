Tutorial
========

Download Template Database
--------------------------------

PANDORA needs a PDB template database to work.
The database can be quickly retrieved by zenodo by running the command-line tool:

`$ pandora-fetch``


Alternatively, the same function can be done in python:


>>> from PANDORA import Database
>>> Database.install_database()


(Advanced) Generate template Database

You can also generate the database from scratch, downloading and parsing the structures directly from IMGT. This will ensure you to have as many templates as possible, as the quickly-retrievable database will not be re-released often. 

The database can be generated with the command-line tool:

`pandora-create``


Or with the pythoncode below: 


>>> ## import requested modules
>>> from PANDORA import Database
>>> 
>>> ## A. Create local Database
>>> db = Database.Database()
>>> db.construct_database(n_jobs=<n_jobs>)

Note 1:  By default, the database generation will use one core only. You can sensitevly speed it up by changing the paramenter --num-cores for pandora-create ot 'n_jobs' for Database.contruct_database().

Note 2: the database is saved by default into `~/PANDORA_databases/default`. It is possible to modify the folder name (`default`) by creating a `config.json` file in the PANDORA installation folder using `data_folder_name` as a key, and the desired folder name as a value, like in the example below:


`{"data_folder_name": "<folder_in_Databases>"}`


p:MHC-I base example
--------------------------------

PANDORA requires at least these information to generate models:

- Peptide sequence
- MHC allele

Steps:

A. Load templates database

B. Create a Target object

C. Generate *n* pMHC models (Default n=20)

Please note that you can specify output directory yourself, otherwise will be generated in the current working directory.

>>> ## import requested modules
>>> from PANDORA import PMHC
>>> from PANDORA import Pandora
>>> from PANDORA import Database
>>>
>>> ## A. Create local Database
>>> db = Database.load()
>>>
>>> ## B. Create Target object
>>> target = PMHC.Target(
>>>     allele_type='HLA-A*02:01'
>>>     peptide='LLFGYPVYV',
>>>     anchors = [2,9])
>>>
>>> ## C. Perform modelling
>>> case = Pandora.Pandora(target, db)
>>> case.model()

Increased loop models
--------------------------------

A user might want to increase the sampling or the sampling quality of the modelling.
To do so, there are two options that can be changed in the Pandora.model function: n_loop_models and loop_refinement.
The first option increases the number of generated loop models, while the second option controls the type of loop refinement MODELLER will perform.

>>> from PANDORA import PMHC
>>> from PANDORA import Pandora
>>> from PANDORA import Database
>>>
>>> ## A. load the pre-generated templates database
>>> db = Database.load()
>>>
>>> ## B. Create Target object
>>> target = PMHC.Target(id='myTestCase'
>>>     allele_type = 'HLA-B*5301',
>>>     peptide = 'TPYDINQML')
>>>
>>> ## C. Perform modelling
>>> case = Pandora.Pandora(target, db)
>>> case.model(n_loop_models=100)  # Generates 100 models

or, to perform a slower and more accurate loop refinement:

>>> case.model(n_loop_models=100)  # Generates 100 models


Run PANDORA on multiple cases (running in parallel on multiple cores)
---------------------------------------------------------------------------

PANDORA can model more than one pMHC complex in parallel. You need to provide the following information in a *.tsv* file:

- *Peptide sequence,  Allele name*


>>> from PANDORA import Pandora
>>> from PANDORA import Database
>>> from PANDORA import Wrapper
>>>
>>> ## A. load the pre-generated templates database
>>> db = Database.load()
>>>
>>> ## B. Create the wrapper object
>>> wrap =  Wrapper(data_file='datafile.tsv', database=db, n_jobs=128)


p:MHC-II base example
-------------------------------------------

>>> from PANDORA import PMHC
>>> from PANDORA import Pandora
>>> from PANDORA import Database
>>>
>>> ## A. Load pre-generated template database
>>> db = Database.load()
>>>
>>> target = PMHC.Target(
>>>     MHC_class='II',
>>>     allele_type = [ 'HLA-DRA*0101', 'HLA-DRB1*0101'],
>>>     peptide = 'GELIGILNAAKVPAD',
>>>     anchors = [4, 7, 9, 12])
>>>
>>> case = Pandora.Pandora(target, db)
>>> case.model()
