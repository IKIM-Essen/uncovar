# trimming of 5' and 3' end
# trimming of bad reads (low qual, too short, too many N)
# cutout adapters
# able to parallelize
rule fastp_pe:
    input:
        # takes both reads per sample
        sample=get_fastqs,
    output:
        # writes out trimmed reads
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
