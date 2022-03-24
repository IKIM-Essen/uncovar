# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
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


# TODO Aggreate assembly_megahit_pe and assembly_megahit_se in one rule.
rule assembly_megahit_pe:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        temp(
            "results/{date}/assembly/{sample}/megahit-{preset}-pe/{sample}.contigs.fasta"
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
        " mv {params.outdir}/final.contigs.fa {output})"
        " > {log} 2>&1"


rule assembly_megahit_se:
    input:
        get_reads_after_qc,
    output:
        temp(
            "results/{date}/assembly/{sample}/megahit-{preset}-se/{sample}.contigs.fasta"
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
        "(megahit -r {input} {params.preset} --out-dir {params.outdir} -f && "
        " mv {params.outdir}/final.contigs.fa {output})"
        " > {log} 2>&1"


# TODO Aggreate assembly_spades_pe and assembly_spades_se in one rule.
rule assembly_spades_pe:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        temp(
            "results/{date}/assembly/{sample}/{spadesflavor}-pe/{sample}.contigs.fasta"
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
        " if [ -f {params.outdir}/raw_contigs.fasta ]; then mv {params.outdir}/raw_contigs.fasta {output}; else mv {params.outdir}/contigs.fasta {output}; fi )"
        " > {log} 2>&1"


rule spades_assemble_se:
    input:
        get_reads_after_qc,
    output:
        temp("results/{date}/assembly/{sample}/spades-se/{sample}.contigs.fasta"),
    log:
        "logs/{date}/spades/se/{sample}.log",
    conda:
        "../envs/spades.yaml"
    params:
        outdir=get_output_dir,
    threads: 8
    shell:
        "if [ -d {params.outdir} ]; then rm -Rf {params.outdir}; fi && "
        "(spades.py --corona -s {input} -o {params.outdir} -t {threads} && "
        " mv {params.outdir}/raw_contigs.fasta {output})"
        " > {log} 2>&1"


rule check_contigs:
    input:
        get_contigs,
    output:
        temp("results/{date}/contigs/checked/{sample}.fasta"),
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
        temp("results/{date}/contigs/ordered/{sample}.fasta"),
    log:
        "logs/{date}/ragoo/{sample}_cleaned.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/ragoo-remove-chr0.py"


# polish illumina de novo assembly
rule assembly_polishing_illumina:
    input:
        fasta="results/{date}/contigs/ordered/{sample}.fasta",
        bcf="results/{date}/filtered-calls/ref~{sample}/{sample}.clonal.nofilter.orf.bcf",
        bcfidx="results/{date}/filtered-calls/ref~{sample}/{sample}.clonal.nofilter.orf.bcf.csi",
    output:
        report(
            "results/{date}/polishing/bcftools-illumina/{sample}.fasta",
            category="4. Sequences",
            subcategory="1. De Novo Assembled Sequences",
            caption="../report/assembly_illumina.rst",
            labels={"sample": "{sample}"},
        ),
    log:
        "logs/{date}/bcftools-consensus-illumina/{sample}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools consensus -f {input.fasta} {input.bcf} > {output} 2> {log}"


# polish ont de novo assembly
rule assembly_polishing_ont:
    input:
        fasta="results/{date}/norm_trim_raw_reads/{sample}/{sample}.cap.clip.fasta",
        reference="results/{date}/contigs/ordered/{sample}.fasta",
    output:
        report(
            "results/{date}/polishing/medaka/{sample}/consensus.fasta",
            category="4. Sequences",
            subcategory="1. De Novo Assembled Sequences",
            caption="../report/assembly_ont.rst",
            labels={"sample": "{sample}"},
        ),
    log:
        "logs/{date}/medaka/consensus/{sample}.log",
    params:
        outdir=get_output_dir,
        model=config["assembly"]["oxford nanopore"]["medaka_model"],
    conda:
        "../envs/medaka.yaml"
    threads: 4
    shell:
        "medaka_consensus -f -i {input.fasta} -o {params.outdir} -d {input.reference} -t {threads} > {log} 2>&1"


rule aggregate_polished_de_novo_sequences:
    input:
        get_polished_sequence,
    output:
        temp("results/{date}/contigs/polished/{sample}.fasta"),
    log:
        "logs/{date}/aggregate_polished_de_novo_sequences/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cp {input} {output} 2> {log}"


rule aggregate_fallback_sequences:
    input:
        get_fallback_sequence,
    output:
        temp("results/{date}/contigs/fallback/{sample}.fasta"),
    log:
        "logs/{date}/aggregate_fallback_sequences/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cp {input} {output} 2> {log}"


rule align_contigs:
    input:
        target="resources/genomes/main.fasta",
        query=get_quast_fastas,
    output:
        temp("results/{date}/aligned/ref~main/{stage}~{sample}.bam"),
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
