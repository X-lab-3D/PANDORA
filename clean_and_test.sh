$clean 
rm -r PANDORA_files/outputs/bm_20210130/*
qsub -q all.q@narrativum.umcn.nl cluster_run_scripts/octarine_run_benchmark.sh bm_20210130 1 1 128
watch --interval=1 qstat -f
