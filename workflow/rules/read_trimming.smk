rule fastp:
    input:
        get_fastqs,
    output:
        fastq1="results/trimmed/{sample}.1.fastq.gz",
        fastq2="results/trimmed/{sample}.2.fastq.gz",
        qc="results/trimmed/{sample}.qc.txt",
    log:
        "logs/fastp/{sample}.log",
    threads: 8
    wrapper:
        "0.70.0/bio/fastp"
