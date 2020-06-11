#!/bin/bash
#SBATCH -J database_download_test
#SBATCH -n 1 -c 80
#SBATCH -t 10

srun -n 1 python ./parallel_neoantigen_model.py test 0 8
