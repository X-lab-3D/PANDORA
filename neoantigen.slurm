#!/usr/bin/bash
#SBATCH -p short
#SBATCH -J job$SLURM_ARRAY_TASK_ID
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH -t 10:00
#SBATCH -o /home/dmarz/pMHC_Modelling/outputs/logs/job$SLURM_ARRAY_TASK_ID.slurm.out
#SBATCH -e /home/dmarz/pMHC_Modelling/outputs/logs/job$SLURM_ARRAY_TASK_ID.slurm.err

    start=`date +%s`

    python ./neoantigen_model.py test $SLURM_ARRAY_TASK_ID $1

    end=`date +%s`

    runtime=$((end-start))
    echo
    echo "total runtime: $runtime sec"

