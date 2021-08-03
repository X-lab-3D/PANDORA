# PANDORA
## (Peptide ANchored mODelling fRAmework)

## Installation

To install PANDORA and all the required dependencies, ```pip3``` is required.
To create all the required directories and install all dependencies:
1. set:
> alias KEY_MODELLER='XXXX'
where 'XXXX' is your MODELLER license key
2. run:
> python install.py


# Dependencies

Python 3

BioPython

Modeller 9.23 or later

pdb_tools

muscle

pdb2sql (only for RMSD calculation)



## PANDORA content

# Contacts

# Database
Contains all the parsing scripts to create database object, retrieve IDs from IMGT, download and clean the structure.

# Pandora
Main class necessary for modelling

# PMHC
Contains PMHC class, superclass of Template, Target and Model 

# run
Running scripts to be removed before the last commit

## PANDORA_files content
Data. Default output directory


## Example run script

```
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

db = Database.Database()
db.construct_database()

target = PMHC.Target('1A1M',
    db.MHCI_data['1A1M'].allele_type,
    db.MHCI_data['1A1M'].peptide,
    M_chain_seq = db.MHCI_data['1A1M'].M_chain_seq,
    anchors = db.MHCI_data['1A1M'].anchors)

mod = Pandora.Pandora(target, db)
mod.model(n_models=20, stdev=0.1, seq_based_templ_selection=True, benchmark=False)
```
