# PANDORA
### Peptide ANchored mODelling fRAmework for peptide-MHC complexes


![PANDORA](https://github.com/DarioMarzella/PANDORA/blob/issue_90/flowchart_pMHCI.png?raw=true)

### Contents

- [Overview](#overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Tutorial](#tutorial)
- [License](./LICENSE)
- [Issues](#issues)

## Overview

PANDORA is anchor restrained modelling pipeline for generating peptide-MHC structures.

PANDORA contains multiple functions to pre-process data and it's able to exploit different crucial domain knowledge provided by the user to guide the modelling.


## Dependencies
PANDORA requires MODELLER, python and some python libraries to be installed. You can install the dependencies via the installation process in PANDORA or install from source.
Installation from source provides you the latest version(recommended).

- [Python](https://www.python.org/) 3
- conda
- pip3
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
Enter the cloned directory and then install all the dependencies!
```
cd PANDORA

python install.py
```
#### 3. (Optional) Install NetMHCpan and/or NetMHCIIpan

PANDORA lets the user if he wills to predict peptide's anchor residues instead of using conventional predefined anchor residues for each allele.
In that case you need to install NetMHCpan (for peptide:MHC class I) and/or NetMHCIIpan (for peptide:MHC class II).
You can install from the [source](https://services.healthtech.dtu.dk/software.php) or simply run:
```
python netMHCpan_install.py
```

## Tutorial


#### Example1 : generating a pMHCI complex with peptide sequence and allele name

#### Example2 : Reproducing a pMHCI complex with known experimental PDB structure

```python
from PANDORA.PMHC import PMHC
from PANDORA.Pandora import Pandora
from PANDORA.Database import Database

## 1. Create local Database
db = Database.Database()
db.construct_database()

## 2. Create Target object
target = PMHC.Target('1A1M',
    db.MHCI_data['1A1M'].allele_type,
    db.MHCI_data['1A1M'].peptide,
    M_chain_seq = db.MHCI_data['1A1M'].M_chain_seq,
    anchors = db.MHCI_data['1A1M'].anchors)

## 3. Perform modelling
mod = Pandora.Pandora(target, db)
mod.model(n_models=20, stdev=0.1, seq_based_templ_selection=True, benchmark=False)
```
## File Structure

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
