#!/bin/bash
#$ -l h_rt=00:10:00
#$ -cwd
#$ -V

## Usage example in array: qsub -q all.q@narrativum.umcn.nl <array of jobs [-t 1-12]> octarine_run_neoantigen.sh <output directory> <n_tasks> <n_cores>
## python ./neoantigen_model.py $1 $SGE_TASK_ID $3 $4 > outputs/logs/log_bm_$1_$SGE_TASK_ID.log

## Usage example without array: qsub -q all.q@narrativum.umcn.nl octarine_run_neoantigen.sh <output directory> <1> <1> <n_cores>
python ./neoantigen_model.py $1 $2 $3 $4 > outputs/logs/log_$1.log
