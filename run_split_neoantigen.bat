#!/bin/bash
#SBATCH -J na_split_test
#SBATCH -n 1 -c 24
#SBATCH -t 5

srun -n 1 python ./parallel_neoantigen_model.py test $SLURM_ARRAY_TASK_ID $argv[1] >log_$SLURM_ARRAY_TASK_ID.log 
