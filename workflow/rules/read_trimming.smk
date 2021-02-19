rule fastp_pe:
    input:
        sample=get_fastqs,
    output:
        trimmed=[
            "results/trimmed/{sample}.1.fastq.gz",
            "results/trimmed/{sample}.2.fastq.gz",
        ],
        html="results/trimmed/{sample}.html",
        json="results/trimmed/{sample}.fastp.json",
    log:
        "logs/fastp/{sample}.log",
    params:
        # adapters=get_adapters,
        adapters=config["adapters"]["illumina"],
        extra="--qualified_quality_phred {} ".format(
            config["RKI-quality-criteria"]["min-PHRED"]
        ) + "--length_required {}".format(
            config["RKI-quality-criteria"]["min-length-reads"]
        ),
    threads: 2
    wrapper:
        "0.70.0/bio/fastp"
