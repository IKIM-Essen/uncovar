rule fastp_pe:
    input:
        sample=get_fastqs,
    output:
        trimmed=[
            "results/trimmed/{sample}.1.fastq.gz",
            "results/trimmed/{sample}.2.fastq.gz",
        ],
        html="results/trimmed/{sample}.html",
        json="results/trimmed/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    params:
        # adapters=get_adapters,
        adapters=config["adapters"]["illumina"],
        extra="",
    threads: 8
    wrapper:
        "0.70.0/bio/fastp"