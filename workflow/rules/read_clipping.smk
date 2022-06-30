# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule samtools_sort:
    input:
        get_samtools_sort_input,
    output:
        temp("results/{date}/read-sorted/{read_type}~{sorted_by}/{sample}.{stage}.bam"),
    params:
        extra=(
            lambda wildcards: "-n -m 4G" if wildcards.sorted_by == "name" else "-m 4G"
        ),
        tmp_dir="/tmp/",
    log:
        "logs/{date}/sort-bam/{read_type}~{sorted_by}/{sample}.{stage}.log",
    threads: 8
    wrapper:
        "0.74.0/bio/samtools/sort"


rule bed_to_tsv:
    input:
        check_bed_for_URL(config["preprocessing"]["amplicon-primers"]),
    output:
        "resources/primer.tsv",
    log:
        "logs/bed-to-bedpe.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/bed-to-tsv.py"


rule trim_primers_fgbio:
    input:
        bam="results/{date}/read-sorted/{read_type}~position/{sample}.initial.bam",
        bai="results/{date}/read-sorted/{read_type}~position/{sample}.initial.bam.bai",
        ref="resources/genomes/{reference}.fasta".format(
            reference=config["preprocessing"]["amplicon-reference"]
        ),
        primers="resources/primer.tsv",
    output:
        "results/{date}/read-clipping/hardclipped/{read_type}/{sample}/{sample}.bam",
    log:
        "logs/{date}/fgbio_primer_clipping/{read_type}/{sample}.log",
    conda:
        "../envs/fgbio.yaml"
    shell:
        "fgbio TrimPrimers -i {input.bam} -p {input.primers} -o {output} -H true -r {input.ref} > {log} 2>&1"


rule filter_bam_fgbio:
    input:
        bam="results/{date}/read-clipping/hardclipped/{read_type}/{sample}/{sample}.bam",
    output:
        "results/{date}/read-clipping/hc_filtered/{read_type}/{sample}/{sample}.bam",
    log:
        "logs/{date}/fgbio_filter_bam/{read_type}/{sample}.log",
    conda:
        "../envs/fgbio.yaml"
    shell:
        "fgbio FilterBam -i {input.bam} -o {output} --min-insert-size 100 --remove-single-end-mappings > {log} 2>&1"


rule samtools_fastq_pe:
    input:
        bam="results/{date}/read-sorted/pe~name/{sample}.hardclipped.bam",
    output:
        fq1=temp("results/{date}/read-clipping/fastq/pe/{sample}.1.fastq.gz"),
        fq2=temp("results/{date}/read-clipping/fastq/pe/{sample}.2.fastq.gz"),
    log:
        "logs/{date}/samtools_fastq/pe/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools fastq -@ {threads} {input.bam} -1 {output.fq1} -2 {output.fq2} 2> {log}"


rule samtools_fastq_se:
    input:
        bam="results/{date}/read-sorted/se~name/{sample}.hardclipped.bam",
    output:
        temp("results/{date}/read-clipping/fastq/se/{sample}.fastq"),
    log:
        "logs/{date}/samtools_fastq/se/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools fastq -@ {threads} {input.bam} > {output} 2> {log}"


rule plot_primer_clipping:
    input:
        unclipped=lambda wildcards: get_input_plotting_primer_clipping(
            wildcards, stage="initial"
        ),
        index_unclipped=lambda wildcards: get_input_plotting_primer_clipping(
            wildcards, stage="initial", suffix=".bai"
        ),
        clipped=lambda wildcards: get_input_plotting_primer_clipping(
            wildcards, stage="hardclipped"
        ),
        index_clipped=lambda wildcards: get_input_plotting_primer_clipping(
            wildcards, stage="hardclipped", suffix=".bai"
        ),
    output:
        plot=report(
            "results/{date}/plots/primer-clipping-intervals.svg",
            caption="../report/amplicon-primer-clipping.rst",
            category="3. Sequencing Details",
            subcategory="4. Correct Amplicon Primer Clipping",
        ),
    params:
        samples=lambda wildcards: get_samples_for_date(wildcards.date),
        bedpe="resources/primer.tsv",
    log:
        "logs/{date}/plot-primer-clipping.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-primer-clipping.py"
