# PANDORA
(Peptide ANchored mODelling fRAmework)

### Contents

- [Overview](#Overview)
- [Installation](#Installation)
- [Quick Tutorial](#Tutorial)
- [License](./LICENSE)
- [Issues](#Issues)

## Overview

![alt text](./image/flowchart_MHCI.png)

PANDORA is a MODELLER-based, anchor restrained modelling pipeline for peptide-MHC structures.

PANDORA contains multiple functions to pre-process data and it's able to exploit different crucial domain knowledge provided by the user to guide its modelling.

## Installation

To install PANDORA and all the required dependencies, ```pip3``` and conda are required.
Also, being based on MODELLER, Pandora needs a functioning MODELLER installation. Please request your MODELLER license at: https://salilab.org/modeller/registration.html
To create all the required directories and install all dependencies:
1. set:
> alias KEY_MODELLER='XXXX'
where 'XXXX' is your MODELLER license key
2. run:
> python install.py


### Dependencies

Python 3
BioPython
Modeller 9.23 or later
pdb_tols
muscle
pdb2sql (only for RMSD calculation)


## Tutorial

### Quick start
```
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

db = Database.Database()
db.construct_database(save = 'a_saved_Pandora_database')

target = PMHC.Target(id = '1A1M',
    allele_type = ['HLA-B*5301', 'HLA-B*5301'],
    peptide = 'TPYDINQML',
    anchors = [2,9])

case = Pandora.Pandora(target, db)
case.model()
```

### Example 1: Create many loop models in a specific directory and calculate the backbone L-RMSD

```
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

db = Database.load('a_saved_Pandora_database')

target = PMHC.Target(id = '1A1M',
    allele_type = ['HLA-B*5301', 'HLA-B*5301'],
    peptide = 'TPYDINQML',
    anchors = [2,9])

case = Pandora.Pandora(target, db)
# to calculate the L-RMSD, the target id must be a PDB id of a template structure present in the database
case.model(output_dir = '/anywhere', n_loop_models = 100, benchmark=True)

```

### Example 2: Model a peptide:MHCI complex with an alpha helix in the peptide

```
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

db = Database.load('a_saved_Pandora_database')

target = PMHC.Target(id = '2YF5',
    allele_type = ['MH1-B*2101', 'MH1-B*2101'],
    peptide = 'TAGQSNYDRL',
    anchors = [2,10],
    helix = ['4', '9'])

case = Pandora.Pandora(target, db)
case.model(helix=target.helix)

```

### Example 3: Predict peptide anchor positions with NetMHCpan

```
from PANDORA.PMHC import PMHC

# Predicting peptide anchor positions requires installation of NetMHCPan (https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1)
target = PMHC.Target(id = '1A1M',
    allele_type = ['HLA-B*5301', 'HLA-B*5301'],
    peptide = 'TPYDINQML)

```

### Example 4: Model many peptide:MHCI complexes in parallel

```
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database
from PANDORA.Wrapper import Wrapper

db = Database.load('a_saved_Pandora_database')

wrap = Pandora.Wrapper()
wrap.create_targets('datafile.tsv', db)
wrap.run_pandora(num_cores = 128)

```

### Example 5: Model a peptide:MHCII complex

```
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

db = Database.load('a_saved_Pandora_database')

target = PMHC.Target(id = '2IAM',
    MHC_class='II',
    allele_type = ['HLA-DRA*0102', 'HLA-DRA*0101', 'HLA-DRB1*0101'],
    peptide = 'GELIGILNAAKVPAD',
    anchors = [4, 7, 9, 12])

case = Pandora.Pandora(target, db)
case.model(stdev=0.2)

```

## Issues

If you have questions or find a bug, please report the issue in the [Github issue channel](https://github.com/DarioMarzella/PANDORA/issues).