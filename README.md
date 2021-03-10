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
import PANDORA

## 1. Create local database
db = PANDORA.Database.construct_database()

## Alternatively, load a pre-built database
dv = PANDORA.Database.load('db.pkl')

## 2. Create target object
target = PANDORA.PMHC.target(PDB_id = '1K5N', MHC_type = 1,
                             allele = ['HLA-B*27:09'], peptide = 'GRFAAAIAK')
                             
## 3. Perform modelling
model = PANDORA.Pandora.model(target, db)
