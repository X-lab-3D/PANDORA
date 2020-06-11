#!/bin/bash
#SBATCH -J na_split_test
#SBATCH -n 1 -c 24
#SBATCH -t 10

srun -n 1 python ./neoantigen_model.py test $SLURM_ARRAY_TASK_ID $1 > outputs/logs/log_$SLURM_ARRAY_TASK_ID.log 
