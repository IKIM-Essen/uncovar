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
    log:
        "logs/{date}/sort-bam/{read_type}~{sorted_by}/{sample}.{stage}.log",
    threads: 8
    wrapper:
        "v1.15.1/bio/samtools/sort"


rule bed_to_bedpe:
    input:
        config["preprocessing"]["amplicon-primers"],
    output:
        "resources/primer.bedpe",
    log:
        "logs/bed-to-bedpe.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/bed-to-bedpe.py"


rule bamclipper:
    input:
        bam="results/{date}/read-sorted/{read_type}~position/{sample}.initial.bam",
        bai="results/{date}/read-sorted/{read_type}~position/{sample}.initial.bam.bai",
        bedpe="resources/primer.bedpe",
    output:
        temp(
            "results/{date}/read-clipping/softclipped/{read_type}/{sample}/{sample}.initial.primerclipped.bam"
        ),
    params:
        output_dir=get_output_dir,
        cwd=lambda w: os.getcwd(),
        bed_path=lambda w, input: os.path.join(os.getcwd(), input.bedpe),
        bam=lambda w, input: os.path.basename(input.bam),
    log:
        "logs/{date}/bamclipper/{read_type}/{sample}.log",
    conda:
        "../envs/bamclipper.yaml"
    threads: 6
    shell:
        "(cp {input.bam} {params.output_dir} &&"
        " cp {input.bai} {params.output_dir} &&"
        " cd {params.output_dir} &&"
        " bamclipper.sh -b {params.bam} -p {params.bed_path} -n {threads} -u 5 -d 5) "
        " > {params.cwd}/{log} 2>&1"


rule fgbio:
    input:
        bam="results/{date}/read-clipping/softclipped/{read_type}/{sample}/{sample}.initial.primerclipped.bam",
        bai="results/{date}/read-clipping/softclipped/{read_type}/{sample}/{sample}.initial.primerclipped.bam.bai",
        ref="resources/genomes/{reference}.fasta".format(
            reference=config["preprocessing"]["amplicon-reference"]
        ),
    output:
        temp(
            "results/{date}/read-clipping/hardclipped/{read_type}/{sample}/{sample}.bam"
        ),
    log:
        "logs/{date}/fgbio/{read_type}/{sample}.log",
    conda:
        "../envs/fgbio.yaml"
    shell:
        "fgbio --sam-validation-stringency=LENIENT ClipBam -i {input.bam} -o {output} -H true -r {input.ref} > {log} 2>&1"


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
        bedpe="resources/primer.bedpe",
    log:
        "logs/{date}/plot-primer-clipping.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-primer-clipping.py"
