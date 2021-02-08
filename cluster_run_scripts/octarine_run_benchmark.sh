#!/bin/bash
#$ -l h_rt=03:00:00
#$ -cwd
#$ -V

## Usage example in array: qsub -q all.q@narrativum.umcn.nl <array of jobs [-t 1-12]> octarine_run_benchmark.sh <output directory> <n_tasks> <n_cores>
## python ./benchmark_model_gap_ali.py $1 $SGE_TASK_ID $3 $4 > outputs/logs/log_bm_$1_$SGE_TASK_ID.log

## Usage example without array: qsub -q all.q@narrativum.umcn.nl octarine_run_benchmark.sh <output directory> <1> <1> <n_cores>
python ./PANDORA/run/benchmark_model.py $1 $2 $3 $4 > ./PANDORA_files/outputs/logs/$1.log
