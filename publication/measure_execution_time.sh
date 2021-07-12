#!/bin/bash
# Script to measure the execution time of the workflow
# Requires an GISAID API token


# number of reads to generate
NO_READS=300000

# number of snakemake executions
NO_ITER=10

# Outfile for time measurements
OUTPUT_LOG=time.log

# Outfile for snakemake log
SNAKE_LOG=snakemake.log

function request_siumlation_samples {
    snakemake --cores all --use-conda results/2021-07-12/tables/strain-genomes.txt
    snakemake --cores all --use-conda resources/benchmarking/B.1.1.7/reads.1.fastq.gz \
        resources/benchmarking/B.1.351/reads.1.fastq.gz 
}

function set_sample_sheet {
    echo "sample_name,fq1,fq2,run_id,is_amplicon_data" > config/pep/samples.csv
    echo "1,resources/benchmarking/B.1.1.7/reads.1.fastq.gz,resources/benchmarking/B.1.1.7/reads.2.fastq.gz,2021-07-12,0" >> config/pep/samples.csv
}

function run_snakemake {
    snakemake --cores all --use-conda
}

# set number of reads to generate
sed -i "s/max_reads: [0-9]*/max_reads: $NO_READS/" config/config.yaml

# simulate samples
echo -e "\nSimulating samples..."
request_siumlation_samples

# update sample sheet accordingly
echo -e "\nUpdating sample sheet..."
set_sample_sheet

#  remove old logs
if [ -f "$OUTPUT_LOG" ] ; then
    rm "$OUTPUT_LOG"
fi

if [ -f "$SNAKE_LOG" ] ; then
    rm "$SNAKE_LOG"
fi

# snakemake time execution loop
echo -e "Starting time measurement"
COUNTER=0
while [  $COUNTER -lt $NO_ITER ]; do
    let COUNTER=COUNTER+1 
    echo "Starting run $COUNTER"
    echo -e "\nRun $COUNTER" >> "$OUTPUT_LOG"

    #  remove results between runs
    if [ -d "results/" ]; then rm -Rf results/; fi

    # time snakemake execution
    { time run_snakemake 2>> $SNAKE_LOG ; } 2>> "$OUTPUT_LOG"
done

