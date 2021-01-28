rule cutadapt:
    input:
        get_fastqs,
    output:
        fastq1="results/trimmed/{sample}.1.fastq.gz",
        fastq2="results/trimmed/{sample}.2.fastq.gz",
        qc="results/trimmed/{sample}.qc.txt",
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="",
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="--minimum-length 1 -q 20",
    log:
        "logs/cutadapt/{sample}.log",
    threads: 8
    wrapper:
        "0.69.0/bio/cutadapt/pe"
