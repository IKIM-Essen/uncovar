# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule samtools_sort:
    input:
        get_samtools_sort_input
    output:
        "results/{date}/read-sorted/{read_type}/{sample}.{stage}.bam",
    params:
        extra="-m 4G",
        tmp_dir="/tmp/",
    log:
        "logs/{date}/sort-bam/{read_type}/{sample}.{stage}.log",
    threads: 8
    wrapper:
        "0.74.0/bio/samtools/sort"


rule bamclipper:
    input:
        bam="results/{date}/read-sorted/{read_type}/{sample}.initial.bam",
        bamidx="results/{date}/read-sorted/{read_type}/{sample}.initial.bam.bai",
        bed=config["adapters"]["amplicon-primers"]
    output:
        temp("results/{date}/read-clipping/softclipped/{read_type}/{sample}/{sample}.initial.primerclipped.bam"),
    params:
        output_dir=get_output_dir,
        cwd=lambda w: os.getcwd(),
        bed_path=lambda w, input: os.path.join(os.getcwd(), input.bed),
        bam=lambda w, input: os.path.basename(input.bam),
    log:
        "logs/{date}/bamclipper/{read_type}/{sample}.log"
    conda:
        "../envs/bamclipper.yaml"
    threads: 6
    shell:
        "(cp {input.bam} {params.output_dir} &&"
        " cp {input.bamidx} {params.output_dir} &&"
        " cd {params.output_dir} &&"
        " bamclipper.sh -b {params.bam} -p {params.bed_path} -n {threads}) "
        " > {params.cwd}/{log} 2>&1"


rule fgbio:
    input:
        bam="results/{date}/read-clipping/softclipped/{read_type}/{sample}/{sample}.initial.primerclipped.bam",
        ref="resources/genomes/{reference}.fasta".format(
            reference=config["adapters"]["amplicon-reference"]
        ),
    output:
        "results/{date}/read-clipping/hardclipped/{read_type}/{sample}/{sample}.bam",
    log:
        "logs/{date}/fgbio/{read_type}/{sample}.log"
    conda:
        "../envs/fgbio.yaml"
    shell:
        "fgbio --sam-validation-stringency=LENIENT ClipBam -i {input.bam} -o {output} -H true -r {input.ref} > {log} 2>&1"


rule samtools_fastq_pe:
    input:
        "results/{date}/read-sorted/{read_type}/{sample}.hardclipped.bam",
    output:
        fq1 = "results/{date}/read-clipping/fastq/se/{sample}.1.fastq.gz",
        fq2 = "results/{date}/read-clipping/fastq/pe/{sample}.2.fastq.gz"
    log:
        "logs/{date}/samtools_fastq/se/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools fastq -@ {threads} {input} -1 {output.fq1} -2 {output.fq2} 2> {log}"


rule samtools_fastq_se:
    input:
        "results/{date}/read-sorted/se/{sample}.hardclipped.bam",
    output:
        "results/{date}/read-clipping/fastq/se/{sample}.fastq"
    log:
        "logs/{date}/samtools_fastq/se/{sample}.log"
    conda:
        "../envs/samtools.yaml"
    threads: 4
    shell:
        "samtools fastq -@ {threads} {input} > {output} 2> {log}"


rule plot_primer_clipping:
    input:
        unclipped=expand_samples_for_date_amplicon(
            "results/{{date}}/read-sorted/pe/{sample}.bam"
        ),
        index_unclipped=expand_samples_for_date_amplicon(
            "results/{{date}}/read-sorted/pe/{sample}.bam.bai"
        ),
        clipped=expand_samples_for_date_amplicon(
            "results/{{date}}/read-sorted/pe/{sample}.primerclipped.hard.sorted.bam"
        ),
        index_clipped=expand_samples_for_date_amplicon(
            "results/{{date}}/read-sorted/pe/{sample}.primerclipped.hard.sorted.bam.bai"
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
        bed=config["adapters"]["amplicon-primers"],
    log:
        "logs/{date}/plot-primer-clipping.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-primer-clipping.py"
