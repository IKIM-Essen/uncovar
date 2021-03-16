rule fastp_pe:
    input:
        sample=get_fastqs,
    output:
        trimmed=[
            "results/{date}/trimmed/{sample}.1.fastq.gz",
            "results/{date}/trimmed/{sample}.2.fastq.gz",
        ],
        html="results/{date}/trimmed/{sample}.html",
        json="results/{date}/trimmed/{sample}.fastp.json",
    log:
        "logs/{date}/fastp/{sample}.log",
    params:
        # adapters=get_adapters,
        adapters=config["adapters"]["illumina"],
        extra="--qualified_quality_phred {} ".format(
            config["RKI-quality-criteria"]["illumina"]["min-PHRED"]
        ) + "--length_required {}".format(
            config["RKI-quality-criteria"]["illumina"]["min-length-reads"]
        ),
    threads: 2
    wrapper:
        "0.70.0/bio/fastp"


rule sum_softclips:
    input:
        "results/{date}/mapped/ref~main/{sample}.bam",
    output:
        "results/{date}/sum-softclips/{sample}.txt",
    log:
        "logs/{date}/sum-softclips/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize-softclips.py"


rule agg_softclips:
    input:
        expand(
            "results/2021-03-13/sum-softclips/{sample}.txt",
            sample=get_samples_for_date("2021-03-13"),
        ),
