#!/bin/bash
#$ -l h_rt=03:00:00
#$ -cwd
#$ -V

## Usage example without array: qsub -q all.q@narrativum.umcn.nl oc_run_rmsds_pdb2sql.py <n_cores>
python ../../tools/add_all_rmsds_pdb2sql.py $1  > ../logs/rmsd_pdb2sql.log

