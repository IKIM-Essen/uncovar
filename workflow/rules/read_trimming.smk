rule fastp_pe:
    input:
        get_fastqs,
    output:
        trimmed=["results/trimmed/{sample}.1.fastq.gz", "results/trimmed/{sample}.2.fastq.gz"],
        html="results/trimmed/{sample}.html",
        json="results/trimmed/{sample}.json"
    log:
        "logs/fastp/{sample}.log"
    params:
        adapters="--adapter_sequence ACGGCTAGCTA --adapter_sequence_r2 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
        extra=""
    threads: 8
    wrapper:
        "0.70.0/bio/fastp"