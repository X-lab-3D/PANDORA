# PANDORA

![Build](https://github.com/X-lab-3D/PANDORA/actions/workflows/main.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/X-lab-3D/PANDORA/badge.svg?branch=master)](https://coveralls.io/github/X-lab-3D/PANDORA?branch=master)

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
- [muscle](https://anaconda.org/bioconda/muscle)
- [Modeller](https://anaconda.org/salilab/modeller) 9.23 or later
- [pdb2sql](https://github.com/DeepRank/pdb2sql) (Optional, only for RMSD calculation)
- [NetMHCpan](https://services.healthtech.dtu.dk/software.php) (Optional, only if user wants to predict peptide:MHC class I anchors)
- [NetMHCIIpan](https://services.healthtech.dtu.dk/software.php) (Optional, only if user wants to predict peptide:MHC class II anchors)

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
git clone https://github.com/X-lab-3D/PANDORA.git
```
Enter the cloned directory and then install the dependencies!
```
cd PANDORA

python install.py
```
#### 3. (Optional) Install NetMHCpan and/or NetMHCIIpan

PANDORA lets the user if he wants to predict peptide's anchor residues instead of using conventional predefined anchor residues.
In that case you need to install NetMHCpan (for peptide:MHC class I) and/or NetMHCIIpan (for peptide:MHC class II).
To install, you can simply run:
```
python netMHCpan_install.py
```

## Tutorial

#### Example 1 : Generating a peptide:MHC complex given the peptide sequence
PANDORA requires at least these information to generate models:
- Peptide sequence
- MHC allele

Steps:

A. The database of all templates need to be generated (retrieving all available pMHC PDBs in [IMGT](http://www.imgt.org/3Dstructure-DB/) database).
   We strongly recommended to save the database once (set argument *save=<your_database_name>*), to skip downloading all templates again for later usage.

B. Creating a Template object based on the given target information

C. Generating *n* number of pMHC models (Default n=20)

Please note that you can specify output directory yourself, otherwise will be generated in a default directory
```python
## import requested modules
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Create local Database
db = Database.Database()
db.construct_database(save='pandora_Database')     

## B. Create Target object
target = PMHC.Target(
    allele_type='HLA-A*0201'
    peptide='LLFGYPVYV',
    anchors = [2,9])

## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model()  
```
#### Example 2 : Create multiple loop models in a your given directory
There are some options provided that you can input them as arguments to the functions.

For instance:
- Generate more models for your modelling case
- Specify the output directory yourself
- Give your target a name
- Predict anchors by NetMHCpan

Please note that, if you do not input *anchors* yourself, it will automatically run NetMHCpan to predict anchors.

```python
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. load the pregenerated Database  of all pMHC PDBs as templates
db = Database.load('pandora_Database')   

## B. Create Target object
target = PMHC.Target(id='myTestCase'
    allele_type = ['HLA-B*5301', 'HLA-B*5301'],
    peptide = 'TPYDINQML')

## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model(n_loop_models=100, output_dir = '/your/directory/')  # Generates 100 models
```

#### Example 3 : Benchmark PANDORA on one modelling case
If you want to evaluate the framework on a target with a known experimental structure:
- Provide the PDB ID for the *Target* class
- Set *benchmark=True* for the modelling
  (calculates L-RMSD to show how far the model is from the near-native structure)

```python
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Load pregenerated database of all pMHC PDBs as templates
db = Database.load('pandora_Database')  

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
If you have some domain knowledge of the peptide conformation, whether it forms secondary structures other than loop (Helix/Beta strand), the framework will consider that while modelling the peptide:
```python
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Load pregenerated database of all pMHC PDBs as templates
db = Database.load('pandora_Database')

## B. Create Target object
target = PMHC.Target(
    allele_type = ['MH1-B*2101', 'MH1-B*2101'],
    peptide = 'TAGQSNYDRL',
    anchors = [2,10],
    helix = ['4', '9'])

## C. Perform modelling
case = Pandora.Pandora(target, db)
case.model(helix=target.helix)
```
#### Example 5: Benchmark PANDORA on multiple cases (running in parallel on multiple cores)
PANDORA can model more than one peptide, in parallel. You need to provide the following peptide information in a *.tsv* file:

- *Peptide sequence,  Allele name, PDB ID* (Optional, only used when reproducing models of known peptide:MHC structures)

The Wrapper class is implemented to run PANDORA in parallel on multiple cores.

```python
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database
from PANDORA.Wrapper import Wrapper

## A. Load pregenerated database of all pMHC PDBs as templates
db = Database.load('pandora_Database')

## B. Create the wrapper object
wrap =  Wrapper()

## C. Create all Target Objects based on peptides in the .tsv file
wrap.create_targets('datafile.tsv', db)

## C. Perform modelling
wrap.run_pandora(num_cores=128)
```
#### Example 6: Generating a peptide:MHC class II complex given the peptide sequence
To model a peptide:MHC class II complex, you only need to specify that in *PMHC.Target()* function: as *MHC_class='II'* (By default it is set to model MHC class I).

```python
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## A. Load pregenerated database of all pMHC PDBs as templates
db = Database.load('pandora_Database')

target = PMHC.Target(
    MHC_class='II',
    allele_type = ['HLA-DRA*0102', 'HLA-DRA*0101', 'HLA-DRB1*0101'],
    peptide = 'GELIGILNAAKVPAD',
    anchors = [4, 7, 9, 12])

case = Pandora.Pandora(target, db)
case.model()
```


## Code Design
PANDORA has been implemented in an Object-Oriented Design(OOD). Resulting in a comprehensible and user-friendly framework.

see [Class Diagram](https://github.com/DarioMarzella/PANDORA/blob/master/images/class_diagram.png?raw=true)

## Output

The following file structure is prepared to store the output files for each case. Each modelling case is given a specific name based on target and template ID.

Please note that the modelling results consisting genretaed models by default are stored in *./PANDORA_files/data/outputs/* directory

- Main outputs: *molpdf_DOPE.tsv, *BL*.pdb, modeller.log(
- Input files prepared for modelling: *contacs_*.list, *.ali*
- *.py* files: MODELLER scripts
- MODELLER by product outputs(Generated during the modelling): *D0*, DL*, *IL*.pdb , , *.ini, *.lrsr, *.rsr, *.sch, ...*

```
PANDORA_files
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
