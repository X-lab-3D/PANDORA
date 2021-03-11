# pMHC_modelling

## Installation

To create all the required directories, run:
>>> python install.py

# Dependencies

Python 3 (Preferred)

BioPython

Dill

Modeller 9.23 or later

pdb_tools

muscle

profit (only for RMSD calculation)

pdb2sql (only for RMSD calculation)


## PANDORA modules

# parsing
Contains all the parsing scripts to retrieve IDs from IMGT, download and clean the structure and parse MODELLER log output

# modelling

Contains the scripts necessary for modelling

# tools

Contains mics scripts, mainly for atom contacts and RMSD calculation

# run

Contains the main run scripts for downloading the dataset and perform large scale modelling.

## PANDORA_files

## Benchmark

## Cluster_run_scripts


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
mod.model(n_models=5, stdev=0.1, seq_based_templ_selection=True, benchmark=True)
```
