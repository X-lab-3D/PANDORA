#!/bin/bash
#$ -l h_rt=1:00:00
#$ -cwd
#$ -V 

## Usage example without array: qsub -q all.q@narrativum.umcn.nl -pe smp 64 octarine_run_db.sh
python ./build_db.py