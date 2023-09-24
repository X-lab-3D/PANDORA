#!/bin/bash

# Process and submit jobs based on a given lookup file.
# The file should be generated from the fasta_generator.py script.

JOB_INTERVAL=5  # interval between job submissions
JOB_BATCH_SIZE=25  # number of jobs after which there's a longer pause
LONG_PAUSE=3600  # 1 hour pause in seconds

submit_job() {
    local count=$1
    local hla=$2
    local len=$3

    local fasta_path="/anchor_analysis/fasta_files_r4/${hla}/anchor_${len}mer_input.fa"
    local output_dir="/anchor_analysis/output_files_r4/${hla}/anchor_${len}mer_output/"

    echo $fasta_path
    mkdir -p $output_dir

    LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q research-hpc -a 'docker(susannakiwala/pvactools:1.5.0b)' -J "anchor_analysis_${count}" \
    -M 8000000 -R 'select[mem>=8000] rusage[mem=8000]' -g '/user/h.xia/sat_analysis_20' \
    -oo "/anchor_analysis/logs_r4/sample_${hla}_${len}.out" -e "/anchor_analysis/logs_r4/sample_${hla}_${len}.err" \
    pvacbind run -e $len --iedb-install-directory /opt/iedb --iedb-retries 50 --binding-threshold 500 \
    --allele-specific-binding-thresholds --keep-tmp-files --n-threads 20 $fasta_path "ANCHOR" $hla \
    MHCflurry MHCnuggetsI NetMHC NetMHCcons NetMHCpan PickPocket SMM SMMPMBEC $output_dir

    echo "Job#${count}, HLA:${hla}, Len:${len}"
    echo "Job#${count}, HLA:${hla}, Len:${len}" >> anchor_jobs_matchlist.txt

    # Long pause every JOB_BATCH_SIZE jobs
    if ((count % JOB_BATCH_SIZE == 0)); then
        sleep $LONG_PAUSE
    fi

    sleep $JOB_INTERVAL
}

main() {
    local count=1
    while IFS=, read -r hla len || [[ -n "$line" ]]
    do
        submit_job $count $hla $len
        count=$((count + 1))
    done < "$1"
}

main "$1"
