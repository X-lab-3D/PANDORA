#!/bin/bash
#SBATCH -J database_download_test
#SBATCH -n 1
#SBATCH -t 01:00:00

srun -n 1 python ./download_template_set.py
