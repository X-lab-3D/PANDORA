#!/bin/bash
#$ -l h_rt=5:00:00
#$ -cwd
#$ -V 

## Usage example without array: qsub -pe smp <n_cores> -q all.q@narrativum.umcn.nl octarine_run_mhcii_bm.sh <n_cores>
python ./MHCII_benchmark.py $1