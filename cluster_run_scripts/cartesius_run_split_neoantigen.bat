#!/bin/bash
#SBATCH -p short
#SBATCH -J wt_split_test
#SBATCH -n 1 -c 24
#SBATCH -t 00:25:00

srun -n 1 python ./neoantigen_model.py $1 $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT > outputs/logs/log_$1_$SLURM_ARRAY_TASK_ID.log 
