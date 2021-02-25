#!/bin/bash
##$ -l h_rt=02:00:00
#$ -cwd
#$ -V

## Usage example without array: qsub -q all.q@narrativum.umcn.nl octarine_run_successrate_hitrate.sh
python ../PANDORA/tools/hitrate_successrate.py 
