# PANDORA

![Build](https://github.com/X-lab-3D/PANDORA/actions/workflows/main.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/X-lab-3D/PANDORA/badge.svg?branch=development)](https://coveralls.io/github/X-lab-3D/PANDORA?branch=development)
[![Anaconda-Server Badge](https://anaconda.org/csb-nijmegen/csb-pandora/badges/version.svg)](https://anaconda.org/csb-nijmegen/csb-pandora)
[![Documentation Status](https://readthedocs.org/projects/csb-pandora/badge/?version=latest)](http://csb-pandora.readthedocs.io/?badge=latest)

### Peptide ANchored mODelling fRAmework for peptide-MHC complexes

![PANDORA](https://github.com/DarioMarzella/PANDORA/blob/master/images/flowchart_pMHCI.png?raw=true)

### Contents

- [Overview](#overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Tutorial](#tutorial)
- [Code Design](#code-design)
- [Output](#output)
- [License](./LICENSE)
- [Issues](#issues)

## Overview

PANDORA is anchor restrained modelling pipeline for generating peptide-MHC structures.

It contains multiple functions to pre-process data and it's able to exploit different crucial domain knowledge provided by the user to guide the modelling.

PANDORA documentation can be found at: https://csb-pandora.readthedocs.io/en/latest/


## Dependencies
PANDORA requires MODELLER, python and some python libraries to be installed.
The following installations are required to start PANDORA installation:

- [Python](https://www.python.org/) 3
- conda
- pip3

The installation process will take care of installing the following dependencies (see [Installation](#installation)), no need to install them yourself.

- [BioPython](https://anaconda.org/conda-forge/biopython)
<!-- - [muscle](https://anaconda.org/bioconda/muscle) -->
<!-- - [Modeller](https://anaconda.org/salilab/modeller) 9.23 or later -->
- [pdb2sql](https://github.com/DeepRank/pdb2sql) (Optional, only for RMSD calculation)
- [NetMHCpan](https://services.healthtech.dtu.dk/software.php) (Optional, only if user wants to predict peptide:MHC class I anchors)
- [NetMHCIIpan](https://services.healthtech.dtu.dk/software.php) (Optional, only if user wants to predict peptide:MHC class II anchors)


## Installation
### Conda Installation (AVAILABLE ONLY FOR MASTER BRANCH)

#### 1. Get a Modeller Key License:
Prior to PANDORA installation, you need to first activate MODELLER's license. Please request MODELLER license at: https://salilab.org/modeller/registration.html

Replace XXXX with your MODELLER License key and run the command:

```
alias KEY_MODELLER='XXXX'
```

#### 2. Install PANDORA

Install with conda:
```
conda install -c csb-nijmegen csb-pandora -c salilab -c bioconda
```

## GitHub / Pypi installation

#### 1. Install Modeller:
Prior to PANDORA installation, you need to first activate MODELLER's license. Please request MODELLER license at: https://salilab.org/modeller/registration.html

Replace XXXX with your MODELLER License key and run the command:
```
alias KEY_MODELLER='XXXX'
```

Then Install MODELLER with:
```
conda install -y -c salilab modeller
```

#### 2. Install Other dependencies
PANDORA relies on muscle (https://anaconda.org/bioconda/muscle) and blast (https://anaconda.org/bioconda/blast) that can be both installed via bioconda.


```
conda install -c bioconda muscle=5.1 blast=2.10
```
For some HPC systems the conda blast installation might not work due to missing libraries. For those cases, you can download blast from their release page and install it (make shure you add it to your PATH, otherwise PANDORA will not be able to find it): https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/

#### 3. Install PANDORA

Clone the repository:
```
git clone https://github.com/X-lab-3D/PANDORA.git
```
Enter the cloned directory, switch to development branch and install the package:
```
cd PANDORA
git checkout development
pip install -e .

```


## Generate / download template Database

PANDORA needs a PDB template database to work. The suggested way is to generate it from scratch, to have as many templates as possible. 
All the structures are retrieved from  [IMGT](http://www.imgt.org/3Dstructure-DB/) database.

You can use the code below to generate your database. By default, the database generation will use all the available cores for the process. You can change it by changing the paramente 'n_jobs' for Database.contruct_database().

```python
## import requested modules
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Create local Database
db = Database.Database()
db.construct_database()
```

Note 1: in the case of generation of template database, data are saved by default into `Databases/default`. It is possible to modify the folder name (`default`) by creating a `config.json` file in the current working directory using `data_folder_name` as a key, and the desired folder name as a value:
```
{"data_folder_name": "<folder_in_Databases>"}
```

Note 2 (For master branch only, conda release v0.9.0): You can download the pre-made database from https://github.com/X-lab-3D/PANDORA_database (pMHC I only, generated on 23/03/2021) and follow the [instructions](https://github.com/X-lab-3D/PANDORA_database/blob/main/README.md). Please be sure you re-path your database as explained in the instructions.

## (Optional) Install NetMHCpan and/or NetMHCIIpan

PANDORA lets the user to predict peptide's anchor residues instead of using conventional predefined anchor residues.
In that case you need to download [NetMHCpan](https://services.healthtech.dtu.dk/cgi-bin/sw_request) (for peptide:MHC class I) and/or [NetMHCIIpan](https://services.healthtech.dtu.dk/cgi-bin/sw_request) (for peptide:MHC class II).
To install, you can simply run:
```
python netMHCpan_install.py
```

## Tutorial

Are you planning to run multiple pMHC models with similar settings?
Then hop directly to example #5 and try our Wrapper module!

#### Example 1 : Generating a peptide:MHC complex given the peptide sequence
PANDORA requires at least these information to generate models:
- Peptide sequence
- MHC allele

Steps:
A. Load the template database (see installation, point 4)

B. Creating a Template object based on the given target information

C. Generating *n* number of pMHC models (Default *n=20*)

Please note that you can specify output directory yourself, otherwise will be generated in a default directory
```python
## import requested modules
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Create local Database
db = Database.load()

## B. Create Target object
target = PMHC.Target(id = 'myTestCase',
    allele_type = 'HLA-A*0201',
    peptide = 'LLFGYPVYV',
    anchors = [2,9])

## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model()
```
#### Example 2: Run PANDORA Wrapper on multiple cases (running in parallel on multiple cores)

PANDORA can model large batches of peptides in parallel. You need to provide the following peptide information in a *.tsv* or *.csv* file:

- *Peptide sequence,  MHC Allele name / MHC chain sequence*
Note: you can also add various information to your file, including anchors for each case, templates, IDs

The Wrapper class will take care of generating PANDORA target objects and parallelize the modelling on the given number of cores:

```python
from PANDORA.Database import Database
from PANDORA.Wrapper import Wrapper

## A. Load pregenerated database of all pMHC PDBs as templates
db = Database.load()

## B. Create the wrapper object
wrap =  Wrapper.Wrapper()

## C. Create all Target Objects based on peptides in the .tsv file
wrap.create_targets('datafile.tsv', db)

## C. Perform modelling
wrap.run_pandora(num_cores=128)
```
#### Example 3: Create multiple loop models in a your given directory
There are some options provided that you can input them as arguments to the functions.

For instance:
- Generate more models for your modelling case
- Specify the output directory yourself
- Give your target a name
- Predict anchors by NetMHCpan

Please note that, if *anchors* is not specified or *use_netmhcpan* is set to *False*, PANDORA will automatically assign canonical anchors (P2 and PΩ). This can be done automatically only for pMHC-I structures.

```python
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. load the pregenerated Database  of all pMHC PDBs as templates
db = Database.load()

## B. Create Target object
target = PMHC.Target(id = 'myTestCase',
    allele_type = ['HLA-B*5301', 'HLA-B*5302'],
    peptide = 'TPYDINQML',
    use_netmhcpan = True)

## C. Perform modelling
case = Pandora.Pandora(target, db, output_dir = '/your/directory/')
case.model(n_loop_models=100)  # Generates 100 models
```

#### Example 4: Model a peptide:MHCI complex with an alpha helix in the peptide

Input domain secondary structure prediction information (Helix/Beta strand):

```python
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Load pregenerated database of all pMHC PDBs as templates
db = Database.load()

## B. Create Target object
target = PMHC.Target(id = 'myMHCITestCase',
    allele_type = 'MH1-B*2101',
    peptide = 'TAGQSNYDRL',
    anchors = [2,10],
    helix = ['4', '9'])

## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model(helix=target.helix)
```
#### Example 5: Generating a peptide:MHC class II complex given the peptide sequence

To model a peptide:MHC class II complex, you only need to specify that in *PMHC.Target()* function: as *MHC_class='II'* (By default it is set to model MHC class I).

```python
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Load pregenerated database of all pMHC PDBs as templates
db = Database.load()

target = PMHC.Target(id='myMHCIITestCase'
    MHC_class = 'II',
    allele_type = ['HLA-DRA*0102', 'HLA-DRA*0101', 'HLA-DRB1*0101'],
    peptide = 'GELIGILNAAKVPAD',
    anchors = [4, 7, 9, 12])

case = Pandora.Pandora(target, db)
case.model()
```
Note: For MHC II, no canonical anchors can be defined. Therefore the user must either install and use NetMHCIIpan or directly input the anchors positions as *anchors* in *PMHC.Target()*

#### Example 6: Benchmark PANDORA on one modelling case

Evaluate the framework on a target with a known experimental structure:
- Provide the PDB ID for the *Target* class
- Set *benchmark=True* for the modelling
  (calculates L-RMSD to show how far the model is from the near-native structure)

```python
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Load pregenerated database of all pMHC PDBs as templates
db = Database.load()

## B. Create Target object
target = PMHC.Target(id='1A1M',
    allele_type=db.MHCI_data['1A1M'].allele_type,
    peptide=db.MHCI_data['1A1M'].peptide,
    anchors = db.MHCI_data['1A1M'].anchors)
    
## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model(benchmark=True)
```
## Code Design
PANDORA has been implemented in an Object-Oriented Design(OOD). Resulting in a comprehensible and user-friendly framework.

see [Class Diagram](https://github.com/DarioMarzella/PANDORA/blob/master/images/class_diagram.png?raw=true)

## Output

The following file structure is prepared to store the output files for each case. Each modelling case is given a specific name based on target and template ID.

Please note that the modelling results consisting genretaed models by default are stored in *./Databases/default/outputs/* directory

- Main outputs: *molpdf_DOPE.tsv, *BL*.pdb, modeller.log(
- Input files prepared for modelling: *contacs_*.list, *.ali*
- *.py* files: MODELLER scripts
- MODELLER by product outputs(Generated during the modelling): *D0*, DL*, *IL*.pdb , , *.ini, *.lrsr, *.rsr, *.sch, ...*

```
Databases
  └── data
     └── outputs                         Default directory to save output
        └── <target_name>_<template_id>  Each user's modelling case is given a specific name

           ├── molpdf_DOPE.tsv           Ranking all models by molpdf and DOPE modeller's scoring functions
           ├── *BL*.pdb                  Final models
           ├── modeller.log              Printing log file generated by MODELLER, describing modelling steps, or any issues arose along modelling

           ├── *.ali                     Alignment file between template(s) and target used for modelling
           ├── contacts_*.list           Contact restraints

           ├── MyLoop.py                 MODELLER script to set loop modelling parameters for the peptide
           ├── cmd_modeller_ini.py       MODELLER script to generate an initial model to extract restraints from
           ├── cmd_modeller.py           MODELLER script to set the main modelling parameters

           ├── *.ini                     Model generated placing the target atoms at the same coordinate as the template's atoms
           ├── *IL*.pdb                  Initial loop model
           └── ...


```

## Issues

If you have questions or find a bug, please report the issue in the [Github issue channel](https://github.com/X-lab-3D/PANDORA/issues).
