#!/bin/bash
#SBATCH -p normal
#SBATCH -J deepdora_test
#SBATCH -n 1 -c 24
#SBATCH -t 01:00:00

## Usage sbatch --array=1-800 cartesius_run_split_deepdora.bat deepdora
srun -n 1 python ./Deepdora_model.py $1 $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_COUNT 24 > outputs/logs/log_$1_$SLURM_ARRAY_TASK_ID.log 
