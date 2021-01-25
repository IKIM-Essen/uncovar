rule targets:
    input:
        fastq1="results/trimmed/{sample}.1.fastq.gz",
        fastq2="results/trimmed/{sample}.2.fastq.gz",
    output:
        contigs="results/megahit/{sample}.fasta",
    log:
        "logs/megahit/{sample}.log",
    threads: 8
    singularity:
        "docker://folker-snakemake"
# need to determine if we want to run entire workflow in one container?
    shell:
        "somecommand {params.i} {params.o}"
