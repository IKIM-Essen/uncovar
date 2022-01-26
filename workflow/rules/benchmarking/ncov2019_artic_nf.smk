# source: https://github.com/connor-lab/ncov2019-artic-nf
rule ncov2019_artic_nf_illumina_data_prep:
    input:
        get_fastqs,
    output:
        d=directory("resources/ref-data/{sample}"),
        fq1=temp("resources/ref-data/{sample}/{sample}_R1.fastq"),
        fq2=temp("resources/ref-data/{sample}/{sample}_R2.fastq"),
    log:
        "logs/ncov2019_artic_nf_illumina_data_prep/{sample}.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "(mkdir -p {output.d} &&"
        " gzip -d {input[0]} -c > {output.fq1} &&"
        " gzip -d {input[1]} -c > {output.fq2})"
        " 2> {log}"


rule ncov2019_artic_nf_illumina:
    input:
        directory="resources/ref-data/{sample}",
    output:
        consensus="results/benchmarking/ncov2019_artic_nf/illumina/{sample}/{sample}.qc.csv",
    log:
        "logs/ncov2019_artic_nf/{sample}.log",
    threads: 8
    params:
        pipeline="connor-lab/ncov2019-artic-nf",
        revision="v1.3.0",
        profile=["conda"],
        flags="--illumina",
        outdir=lambda w: f"results/benchmarking/ncov2019_artic_nf/illumina/{w.sample}",
        prefix=lambda w: w.sample,
    handover: True
    conda:
        "../../envs/nextflow.yaml"
    script:
        "../../scripts/benchmarking/nextflow.py"


use rule ncov2019_artic_nf_illumina as ncov2019_artic_nf_nanopore_nanopolish with:
    input:
        basecalled_fastq=lambda wildcards: get_fastq_pass_path(wildcards),
        fast5_pass=lambda wildcards: get_fast5_pass_path(wildcards),
        sequencing_summary=lambda wildcards: get_seq_summary(wildcards),
    output:
        consensus="results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{sample}/{sample}.qc.csv",
    log:
        "logs/ncov2019_artic_nf/nanopore/nanopolish/{sample}.log",
    params:
        pipeline="connor-lab/ncov2019-artic-nf",
        revision="v1.3.0",
        profile=["conda"],
        flags="--nanopolish",
        outdir=lambda w: f"results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{w.sample}",
        prefix=lambda w: w.sample,


use rule ncov2019_artic_nf_illumina as bncov2019_artic_nf_nanopore_medaka with:
    input:
        basecalled_fastq=lambda wildcards: get_fastq_pass_path(wildcards),
        fast5_pass=lambda wildcards: get_fast5_pass_path(wildcards),
        sequencing_summary=lambda wildcards: get_seq_summary(wildcards),
    output:
        consensus="results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{sample}/{sample}.qc.csv",
    log:
        "logs/ncov2019_artic_nf/nanopore/medaka/{sample}.log",
    params:
        pipeline="connor-lab/ncov2019-artic-nf",
        revision="v1.3.0",
        profile=["conda"],
        flags="--medaka",
        outdir=lambda w: f"results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{w.sample}",
        prefix=lambda w: w.sample,
