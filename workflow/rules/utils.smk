# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule tabix_index:
    input:
        "{prefix}.{fmt}.gz",
    output:
        temp("{prefix}.{fmt}.gz.tbi"),
    params:
        "-p {fmt}",
    log:
        "logs/tabix-{fmt}/{prefix}.log",
    wrapper:
        "0.70.0/bio/tabix"


rule bam_index:
    input:
        "{prefix}.bam",
    output:
        temp("{prefix}.bam.bai"),
    log:
        "logs/bam-index/{prefix}.log",
    wrapper:
        "0.70.0/bio/samtools/index"


rule bcf_index:
    input:
        "{prefix}.bcf",
    output:
        temp("{prefix}.bcf.csi"),
    log:
        "logs/bcf-index/{prefix}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index {input} 2> {log}"


rule faidx:
    input:
        "{prefix}.fasta",
    output:
        temp("{prefix}.fasta.fai"),
    log:
        "logs/faidx/{prefix}.log",
    wrapper:
        "0.70.0/bio/samtools/faidx"


rule gzip:
    input:
        "{prefix}.fastq",
    output:
        temp("{prefix}.fastq.gz"),
    log:
        "logs/gzip/{prefix}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip --keep {input}"
