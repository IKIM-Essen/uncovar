# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule tabix_index:
    input:
        "{prefix}.{fmt}.gz",
    output:
        "{prefix}.{fmt}.gz.tbi",
    params:
        "-p {fmt}",
    log:
        "logs/tabix-{fmt}/{prefix}.log",
    wrapper:
        "v1.15.1/bio/tabix/index"


rule bam_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/bam-index/{prefix}.log",
    wrapper:
        "v1.15.1/bio/samtools/index"


rule bcf_index:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi",
    log:
        "logs/bcf-index/{prefix}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index {input} 2> {log}"


rule bcf_sort:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.sorted.bcf",
    log:
        "logs/bcf-sort/{prefix}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools sort -O b {input} -o {output} 2> {log}"


rule faidx:
    input:
        "{prefix}.fasta",
    output:
        "{prefix}.fasta.fai",
    log:
        "logs/faidx/{prefix}.log",
    wrapper:
        "v1.15.1/bio/samtools/faidx"


rule gzip:
    input:
        "{prefix}.fastq",
    output:
        "{prefix}.fastq.gz",
    log:
        "logs/gzip/{prefix}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "gzip --keep {input}"
