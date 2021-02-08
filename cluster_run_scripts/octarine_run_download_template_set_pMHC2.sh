#!/bin/bash
##$ -l h_rt=03:00:00
#$ -cwd
#$ -V

## Usage example without array: qsub -q all.q@narrativum.umcn.nl octarine_run_download_template_set.sh
python ./PANDORA/run/download_template_set_pMHC2.py > ./PANDORA_files/outputs/logs/log_test_download.log
