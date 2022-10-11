# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule bwa_index:
    input:
        get_reference(),
    output:
        idx=multiext(
            "results/{date}/bwa/index/ref~{reference}.fasta",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    params:
        prefix=lambda w, output: get_bwa_index_prefix(output),
    log:
        "logs/{date}/bwa-index/ref~{reference}.log",
    wrapper:
        "v1.12.2/bio/bwa/index"


rule bwa_large_index:
    input:
        get_reference(),
    output:
        idx=multiext(
            "resources/bwa/index/ref~{reference}.fasta",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    params:
        prefix=lambda w, output: get_bwa_index_prefix(output),
    log:
        "logs/bwa-index/ref~{reference}.log",
    wrapper:
        "v1.12.2/bio/bwa/index"


rule map_reads:
    input:
        reads=get_reads,
        idx=get_bwa_index,
    output:
        temp("results/{date}/mapped/ref~{reference}/{sample}.bam"),
    log:
        "logs/{date}/bwa-mem/ref~{reference}/{sample}.log",
    params:
        index=lambda w, input: get_bwa_index_prefix(input.idx),
        extra="",
        sorting="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "v1.12.2/bio/bwa/mem"


rule mark_duplicates:
    input:
        bams="results/{date}/mapped/ref~{reference}/{sample}.bam",
    output:
        bam=temp("results/{date}/dedup/ref~{reference}/{sample}.bam"),
        metrics="results/{date}/qc/dedup/ref~{reference}/{sample}.metrics.txt",
    log:
        "logs/{date}/picard/dedup/ref~{reference}/{sample}.log",
    params:
        "",
    wrapper:
        "v1.12.2/bio/picard/markduplicates"


rule samtools_calmd:
    input:
        aln=get_recal_input,
        ref=get_reference(),
    output:
        "results/{date}/recal/ref~{reference}/{sample}.bam",
    log:
        "logs/{date}/samtools-calmd/ref~{reference}/{sample}.log",
    params:
        "-A",
    threads: 8
    wrapper:
        "v1.12.2/bio/samtools/calmd"
