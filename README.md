# PANDORA
### Peptide ANchored mODelling fRAmework for peptide-MHC complexes


![PANDORA](https://github.com/DarioMarzella/PANDORA/blob/issue_90/flowchart_pMHCI.png?raw=true)

### Contents

- [Overview](#overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Tutorial](#tutorial)
- [Code Design](#diagram)
- [Output File structure](#output)
- [License](./LICENSE)
- [Issues](#issues)

## Overview

PANDORA is anchor restrained modelling pipeline for generating peptide-MHC structures.

PANDORA contains multiple functions to pre-process data and it's able to exploit different crucial domain knowledge provided by the user to guide the modelling.


## Dependencies
PANDORA requires MODELLER, python and some python libraries to be installed. 
The following installations are required to start PANDORA installation:

- [Python](https://www.python.org/) 3
- conda
- pip3

PANDORA installation will take care of installing the following dependencies (recommended), no need to install them yourself.

- [BioPython](https://anaconda.org/conda-forge/biopython)
- [muscle](https://anaconda.org/bioconda/muscle)
- [Modeller](https://anaconda.org/salilab/modeller) 9.23 or later
- [pdb_tools](https://github.com/haddocking/pdb-tools)
- [pdb2sql](https://github.com/DeepRank/pdb2sql) (Optional, only for RMSD calculation)
- [NetMHCpan](https://services.healthtech.dtu.dk/software.php) (Optional, only if user wills to predict peptide:MHC class I anchors)
- [NetMHCIIpan](https://services.healthtech.dtu.dk/software.php) (Optional, only if user wills to predict peptide:MHC class II anchors)

## Installation

#### 1. Setup MODELLER License:
Prior to PANDORA installation, you need to first activate MODELLER's license. Please request MODELLER license at: https://salilab.org/modeller/registration.html

Replace XXXX with your MODELLER License key and run the command:
```
alias KEY_MODELLER='XXXX'
```

#### 2. Install PANDORA

Clone the repository:
```
git clone https://github.com/DarioMarzella/PANDORA.git
```
Enter the cloned directory and then install the dependencies!
```
cd PANDORA

python install.py
```
#### 3. (Optional) Install NetMHCpan and/or NetMHCIIpan

PANDORA lets the user if he wills to predict peptide's anchor residues instead of using conventional predefined anchor residues for each allele.
In that case you need to install NetMHCpan (for peptide:MHC class I) and/or NetMHCIIpan (for peptide:MHC class II).
To install, you can simply run:
```
python netMHCpan_install.py
```

## Tutorial


#### Example 1 : Generating a peptide:MHC complex given the peptide sequence
PANDORA requires these information to generate models:
- Peptide sequence
- MHC allele

Steps:

A. The database of all templates need to be generated (retrieving pMHC PDBs from IMGT). 
   Better save the database with your given name, to skip the downloading phase for later usage (recommended).
   
B. Creating a Template object based on the given target information

C. Generating *n* number of pMHC models (Default n=20)


```
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Create local Database
db = Database.Database()
db.construct_database(save='pandora.Database')     

## B. Create Target object
target = PMHC.Target(
    dallele_type='HLA-A*0201'
    peptide='LLFGYPVYV',
    MHC_class='I'
    anchors = [2,9])

## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model()  
```
#### Example 2 : Create many loop models in a specific directory 
There are some options that you can set youself. You can input these as arguments to the functions. For instance if you want to generate 100 models for your modelling case, and also specify the output directory yourself:
```
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Create local Database
db = Database.Database()
db.construct_database(save='pandora.Database')     

## B. Create Target object
target = PMHC.Target(
    allele_type = ['HLA-B*5301', 'HLA-B*5301'],
    peptide = 'TPYDINQML',
    anchors = [2,9])

## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model(n_loop_models=100, output_dir = '/anywhere')  # Generates 100 models
```

#### Example 3 : Reproducing a pMHCI complex with known experimental PDB structure

If the target experimental structure is available, you can provide the PDB ID and set *benchmark=True* to calculate L-RMSD value.

```
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Load pregenerated database of all pMHC PDBs
db = Database.load('pandora.Database')  

## B. Create Target object
target = PMHC.Target('1A1M',
    db.MHCI_data['1A1M'].allele_type,
    db.MHCI_data['1A1M'].peptide,
    anchors = db.MHCI_data['1A1M'].anchors)

## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model(benchmark=True)  
```
#### Example 4: Model a peptide:MHCI complex with an alpha helix in the peptide
If you have some knowledge of the peptide conformation, that it forms secondary structures other than loop, like Helix, or Beta sheet, PANDORA provides you to input these information:
```
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

db = Database.load('a_saved_Pandora_database')

target = PMHC.Target(
    allele_type = ['MH1-B*2101', 'MH1-B*2101'],
    peptide = 'TAGQSNYDRL',
    anchors = [2,10],
    helix = ['4', '9'])

case = Pandora.Pandora(target, db)
case.model(helix=target.helix)
```
#### Example 5: Modelling of many peptide cases
PANDORA can modell more than one peptide, in parallel. You need to provide the following petide information in a csv file, including:

- Peptide sequence,  Allele name, PDB ID (Optional, only used when reproducing models of known peptide:MHC structures)

```
from PANDORA.Wrapper import Wrapper

## A. Load pregenerated database of all pMHC PDBs
db = Database.load('pandora.Database')

## B. Create the wrapper object
wrap =  Wrapper()

## C. Create all Target Objects
wrap.create_targets('datafile.tsv', db, MHC_class='II')

## C. Perform modelling
wrap.run_pandora(num_cores=128)
```
## Diagram
PANDORA is designed in an Object-oriented programming. This provides a comprehensible and user-friendly framework.

Pandora class: represents a user defined modelling case(s)

[Diagram](https://github.com/DarioMarzella/PANDORA/blob/issue_90/class_diagram.png?raw=true)

## Output

The following file structure is prepared to store the Database, PDB files and output data.

Please note that the modelling results consisting genretaed models are stored in *./PANDORA_files/data/outputs/* directory

```
PANDORA_files
  └── data
     ├── csv_pkl_files            Generated Database and csv files containing peptide data to be modelled
     │   └── mhcseqs
     ├── outputs                  Directory to save output (Modelled cases, restraints, alignment file, log file, etc.)
     └── PDBs                     PDB files downloaded from IMGT
           ├── Bad                Problematic PDB files deleted from the databse
           ├── IMGT_retrieved     
           ├── pMHCI              
           └── pMHCII             
```

## Issues

If you have questions or find a bug, please report the issue in the [Github issue channel](https://github.com/DarioMarzella/PANDORA/issues).
