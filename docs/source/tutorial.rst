Tutorial
========


p:MHC-I basic example
--------------------------------

PANDORA requires at least these information to generate models:
- Peptide sequence
- MHC allele

Steps:

A. The database of all templates need to be generated (retrieving all available pMHC PDBs in [IMGT](http://www.imgt.org/3Dstructure-DB/) database).
   We strongly recommended to save the database once (set argument *save=<your_database_name>*), to skip downloading all templates again for later usage.

B. Creating a Template object based on the given target information

C. Generating *n* number of pMHC models (Default n=20)

Please note that you can specify output directory yourself, otherwise will be generated in a default directory

>>> >>> ## import requested modules
>>> from PANDORA.PMHC import PMHC
>>> from PANDORA.Pandora import Pandora
>>> from PANDORA.Database import Database
>>>
>>> ## A. Create local Database
>>> db = Database.Database()
>>> db.construct_database(save='pandora_Database')
>>>
>>> ## B. Create Target object
>>> target = PMHC.Target(
>>>     allele_ty:pe='HLA-A*0201'
>>>     peptide='LLFGYPVYV',
>>>     anchors = [2,9])
>>>
>>> ## C. Perform modelling
>>> case = Pandora.Pandora(target, db)
>>> case.model()

Increased loop models
--------------------------------

There are some options provided that you can input them as arguments to the functions.

For instance:
- Generate more models for your modelling case
- Specify the output directory yourself
- Give your target a name
- Predict anchors by NetMHCpan

Please note that, if you do not input *anchors* yourself, it will automatically run NetMHCpan to predict anchors.


>>> from PANDORA.PMHC import PMHC
>>> from PANDORA.Pandora import Pandora
>>> from PANDORA.Database import Database
>>>
>>> ## A. load the pregenerated Database  of all pMHC PDBs as templates
>>> db = Database.load('pandora_Database')
>>>
>>> ## B. Create Target object
>>> target = PMHC.Target(id='myTestCase'
>>>     allele_type = ['HLA-B*5301', 'HLA-B*5301'],
>>>     peptide = 'TPYDINQML')
>>>
>>> ## C. Perform modelling
>>> case = Pandora.Pandora(target, db)
>>> case.model(n_loop_models=100, output_dir = '/your/directory/')  # Generates 100 models


Benchmark PANDORA on one p:MHC-I case
-------------------------------------

If you want to evaluate the framework on a target with a known experimental structure:
- Provide the PDB ID for the *Target* class
- Set *benchmark=True* for the modelling
(calculates L-RMSD to show how far the model is from the near-native structure)

>>> from PANDORA.PMHC import PMHC
>>> from PANDORA.Pandora import Pandora
>>> from PANDORA.Database import Database
>>>
>>> ## A. Load pregenerated database of all pMHC PDBs as templates
>>> db = Database.load('pandora_Database')
>>>
>>> ## B. Create Target object
>>> target = PMHC.Target('1A1M',
>>>     db.MHCI_data['1A1M'].allele_type,
>>>     db.MHCI_data['1A1M'].peptide,
>>>     anchors = db.MHCI_data['1A1M'].anchors)
>>>
>>> ## C. Perform modelling
>>> case = Pandora.Pandora(target, db)
>>> case.model(benchmark=True)

p:MHC-I complex with an alpha helix in the peptide
--------------------------------------------------

If you have some domain knowledge of the peptide conformation, whether it forms secondary structures other than loop (Helix/Beta strand), the framework will consider that while modelling the peptide:


>>> from PANDORA.PMHC import PMHC
>>> from PANDORA.Pandora import Pandora
>>> from PANDORA.Database import Database
>>>
>>> ## A. Load pregenerated database of all pMHC PDBs as templates
>>> db = Database.load('pandora_Database')
>>>
>>> ## B. Create Target object
>>> target = PMHC.Target(
>>>     allele_type = ['MH1-B*2101', 'MH1-B*2101'],
>>>     peptide = 'TAGQSNYDRL',
>>>     anchors = [2,10],
>>>     helix = ['4', '9'])
>>>
>>> ## C. Perform modelling
>>> case = Pandora.Pandora(target, db)
>>> case.model(helix=target.helix)


Benchmark PANDORA on multiple cases (running in parallel on multiple cores)
---------------------------------------------------------------------------

PANDORA can model more than one peptide, in parallel. You need to provide the following peptide information in a *.tsv* file:

- *Peptide sequence,  Allele name, PDB ID* (Optional, only used when reproducing models of known peptide:MHC structures)

The Wrapper class is implemented to run PANDORA in parallel on multiple cores.


>>> from PANDORA.Pandora import Pandora
>>> from PANDORA.Database import Database
>>> from PANDORA.Wrapper import Wrapper
>>>
>>> ## A. Load pregenerated database of all pMHC PDBs as templates
>>> db = Database.load('pandora_Database')
>>>
>>> ## B. Create the wrapper object
>>> wrap =  Wrapper()
>>>
>>> ## C. Create all Target Objects based on peptides in the .tsv file
>>> wrap.create_targets('datafile.tsv', db)
>>>
>>> ## C. Perform modelling
>>> wrap.run_pandora(num_cores=128)


p:MHC-II complex given the peptide sequence
-------------------------------------------

To model a peptide:MHC class II complex, you only need to specify that in *PMHC.Target()* function: as *MHC_class='II'* (By default it is set to model MHC class I).


>>> from PANDORA.PMHC import PMHC
>>> from PANDORA.Pandora import Pandora
>>> from PANDORA.Database import Database
>>>
>>> ## A. Load pregenerated database of all pMHC PDBs as templates
>>> db = Database.load('pandora_Database')
>>>
>>> target = PMHC.Target(
>>>     MHC_class='II',
>>>     allele_type = ['HLA-DRA*0102', 'HLA-DRA*0101', 'HLA-DRB1*0101'],
>>>     peptide = 'GELIGILNAAKVPAD',
>>>     anchors = [4, 7, 9, 12])
>>>
>>> case = Pandora.Pandora(target, db)
>>> case.model()
