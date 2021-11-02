# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule count_assembly_reads:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
    output:
        read_count=temp("results/{date}/tables/read_pair_counts/{sample}.txt"),
    log:
        "logs/{date}/read_pair_counts/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        '(echo "$(( $(zcat {input.fastq1} | wc -l) / 4 ))" > {output.read_count}) 2> {log}'


rule assembly_megahit:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs=temp(
            "results/{date}/assembly/{sample}/megahit-{preset}/{sample}.contigs.fasta"
        ),
    wildcard_constraints:
        preset="std|meta-large|meta-sensitive",
    log:
        "logs/{date}/megahit-{preset}/{sample}.log",
    params:
        outdir=get_output_dir,
        preset=get_megahit_preset,
    threads: 8
    conda:
        "../envs/megahit.yaml"
    shell:
        "(megahit -1 {input.fastq1} -2 {input.fastq2} {params.preset} --out-dir {params.outdir} -f && "
        " mv {params.outdir}/final.contigs.fa {output.contigs} )"
        " > {log} 2>&1"


rule assembly_spades_pe:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs=temp(
            "results/{date}/assembly/{sample}/{spadesflavor}/{sample}.contigs.fasta"
        ),
    wildcard_constraints:
        spadesflavor="spades|rnaviralspades|metaspades|coronaspades",
    params:
        outdir=get_output_dir,
    log:
        "logs/{date}/{spadesflavor}/pe/{sample}.log",
    conda:
        "../envs/spades.yaml"
    threads: 8
    shell:
        "({wildcards.spadesflavor}.py -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir} -t {threads} && "
        " mv {params.outdir}/*ontigs.fasta {output.contigs})"
        " > {log} 2>&1"


rule spades_assemble_se:
    input:
        "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
    output:
        "results/{date}/assembly/{sample}/spades_se/{sample}.contigs.fasta",
    log:
        "logs/{date}/spades/se/{sample}.log",
    conda:
        "../envs/spades.yaml"
    params:
        outdir=get_output_dir,
    threads: 8
    shell:
        "(spades.py --corona -s {input} -o {params.outdir} -t {threads} && "
        " mv {params.outdir}/raw_contigs.fasta {output} )"
        " > {log} 2>&1"


rule check_contigs:
    input:
        get_contigs,
    output:
        "results/{date}/contigs/checked/{sample}.fasta",
    log:
        "logs/{date}/check_contigs/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/check_contigs.py"


rule order_contigs:
    input:
        contigs="results/{date}/contigs/checked/{sample}.fasta",
        reference="resources/genomes/main.fasta",
    output:
        temp("results/{date}/contigs/ordered-unfiltered/{sample}.fasta"),
    log:
        "logs/{date}/ragoo/{sample}.log",
    params:
        outdir=get_output_dir,
    conda:
        "../envs/ragoo.yaml"
    shadow:
        "minimal"
    shell:
        "(mkdir -p {params.outdir}/{wildcards.sample} && cd {params.outdir}/{wildcards.sample} &&"
        " ragoo.py ../../../../../{input.contigs} ../../../../../{input.reference} &&"
        " cd ../../../../../ && mv {params.outdir}/{wildcards.sample}/ragoo_output/ragoo.fasta {output})"
        " > {log} 2>&1"


rule filter_chr0:
    input:
        "results/{date}/contigs/ordered-unfiltered/{sample}.fasta",
    output:
        "results/{date}/contigs/ordered/{sample}.fasta",
    log:
        "logs/{date}/ragoo/{sample}_cleaned.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/ragoo-remove-chr0.py"


# polish illumina
rule bcftools_consensus:
    input:
        fasta="results/{date}/contigs/ordered/{sample}.fasta",
        bcf="results/{date}/filtered-calls/ref~{sample}/{sample}.clonal.nofilter.bcf",
        bcfidx="results/{date}/filtered-calls/ref~{sample}/{sample}.clonal.nofilter.bcf.csi",
    output:
        report(
            "results/{date}/polishing/bcftools/{sample}.fasta",
            category="4. Assembly",
            caption="../report/assembly_illumina.rst",
        ),
    log:
        "logs/{date}/bcftools-consensus/{sample}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools consensus -f {input.fasta} {input.bcf} > {output} 2> {log}"


# polish ont
rule medaka:
    input:
        fasta="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
        reference="results/{date}/contigs/ordered/{sample}.fasta",
    output:
        report(
            "results/{date}/polishing/medaka/{sample}/consensus.fasta",
            category="4. Assembly",
            caption="../report/assembly_ont.rst",
        ),
    log:
        "logs/{date}/medaka/{sample}.log",
    params:
        outdir=get_output_dir,
        model=config["assembly"]["medaka_model"],
    conda:
        "../envs/medaka.yaml"
    threads: 4
    shell:
        "medaka_consensus -v -i {input.fasta} -o {params.outdir} -d {input.reference} -t {threads} > {log} 2>&1"


rule aggregate_polished_sequences:
    input:
        get_polished_sequence,
    output:
        "results/{date}/contigs/polished/{sample}.fasta",
    log:
        "logs/{date}/aggregate_polished_sequences/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cp {input} {output} 2> {log}"


rule align_contigs:
    input:
        target="resources/genomes/main.fasta",
        query=get_quast_fastas,
    output:
        "results/{date}/aligned/ref~main/{stage}~{sample}.bam",
    log:
        "results/{date}/aligned/ref~main/{stage}~{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax asm5 {input.target} {input.query} -o {output} 2> {log}"


rule quast:
    input:
        fasta=get_quast_fastas,
        bam="results/{date}/aligned/ref~main/{stage}~{sample}.bam",
        reference="resources/genomes/main.fasta",
    output:
        "results/{date}/quast/{stage}/{sample}/report.tsv",
    params:
        outdir=get_output_dir,
    log:
        "logs/{date}/quast/{stage}/{sample}.log",
    conda:
        "../envs/quast.yaml"
    threads: 8
    shell:
        "quast.py --min-contig 1 --threads {threads} -o {params.outdir} -r {input.reference} "
        "--bam {input.bam} {input.fasta} "
        "> {log} 2>&1"
