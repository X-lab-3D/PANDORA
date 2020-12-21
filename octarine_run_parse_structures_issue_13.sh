#!/bin/bash
##$ -l h_rt=03:00:00
#$ -cwd
#$ -V

## Usage example without array: qsub -q all.q@narrativum.umcn.nl octarine_run_parse_structures_issue_13.sh
python ./parse_structures_issue_13.py > outputs/logs/log_test_download.log
