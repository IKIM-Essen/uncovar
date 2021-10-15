# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule fastqc:
    input:
        get_fastqs,
    output:
        html="results/{date}/qc/fastqc/{sample}.html",
        zip="results/{date}/qc/fastqc/{sample}_fastqc.zip",
    log:
        "logs/{date}/fastqc/{sample}.log",
    wrapper:
        "0.69.0/bio/fastqc"


rule multiqc:
    input:
        expand_samples_for_date(
            [
                "results/{{date}}/qc/fastqc/{sample}_fastqc.zip",
                "results/{{date}}/species-diversity-nonhuman/{sample}/{sample}.cleaned.kreport2",
                "results/{{date}}/quast/unpolished/{sample}/report.tsv",
                "results/{{date}}/quast/polished/{sample}/report.tsv",
                "results/{{date}}/qc/samtools_flagstat/{sample}.bam.flagstat",
                "results/{{date}}/qc/dedup/ref~main/{sample}.metrics.txt",
            ]
        ),
        expand_samples_for_date("logs/{{date}}/kallisto_quant/{sample}.log"),
        get_fastp_results,
        get_kraken_output,
    output:
        "results/{date}/qc/multiqc.html",
    params:
        "--config config/multiqc_config.yaml",
        "--title 'Results for data from {date}'",  # Optional: extra parameters for multiqc.
    log:
        "logs/{date}/multiqc.log",
    wrapper:
        "0.69.0/bio/multiqc"


rule multiqc_lab:
    input:
        expand_samples_for_date(
            [
                "results/{{date}}/qc/fastqc/{sample}_fastqc.zip",
                "results/{{date}}/quast/unpolished/{sample}/report.tsv",
            ]
        ),
        get_fastp_results,
        get_kraken_output,
    output:
        report(
            "results/{date}/qc/laboratory/multiqc.html",
            htmlindex="multiqc.html",
            caption="../report/multi-qc-lab.rst",
            category="3. Sequencing Details",
            subcategory="1. Quality Control",
        ),
    params:
        "--config config/multiqc_config_lab.yaml",
        "--title 'Results for data from {date}'",  # Optional: extra parameters for multiqc.
    log:
        "logs/{date}/multiqc.log",
    wrapper:
        "0.69.0/bio/multiqc"


rule samtools_flagstat:
    input:
        "results/{date}/recal/ref~main/{sample}.bam",
    output:
        "results/{date}/qc/samtools_flagstat/{sample}.bam.flagstat",
    log:
        "logs/{date}/samtools/{sample}_flagstat.log",
    wrapper:
        "0.70.0/bio/samtools/flagstat"


rule samtools_depth:
    input:
        get_depth_input,
    output:
        "results/{date}/qc/samtools_depth/{sample}.txt",
    log:
        "logs/{date}/samtools/{sample}_depth.txt",
    conda:
        "../envs/samtools.yaml"
    params:
        ref=config["adapters"]["amplicon-reference"],
    shell:
        "(samtools depth -aH -o {output} {input} && "
        " sed -i 's/{params.ref}.3/{wildcards.sample}/' {output})"
        " 2> {log}"


# analysis of species diversity present BEFORE removing human contamination
rule species_diversity_before_pe:
    input:
        db="resources/minikraken-8GB",
        reads=get_trimmed_reads,
    output:
        classified_reads=temp(
            expand(
                "results/{{date}}/species-diversity/pe/{{sample}}/{{sample}}_{read}.classified.fasta",
                read=[1, 2],
            )
        ),
        unclassified_reads=temp(
            expand(
                "results/{{date}}/species-diversity/pe/{{sample}}/{{sample}}_{read}.unclassified.fasta",
                read=[1, 2],
            )
        ),
        kraken_output=temp(
            "results/{date}/species-diversity/pe/{sample}/{sample}.kraken"
        ),
        report="results/{date}/species-diversity/pe/{sample}/{sample}.uncleaned.kreport2",
    log:
        "logs/{date}/kraken/{sample}.log",
    params:
        classified=lambda w, output: " --unclassified-out "
        + "#".join(output.classified_reads[0].rsplit("_1", 1))
        if is_illumina(w)
        else "",
        unclassified=lambda w, output: "--classified-out "
        + "#".join(output.unclassified_reads[0].rsplit("_1", 1))
        if is_illumina(w)
        else "",
        paired=lambda wildcards: "--paired " if is_illumina(wildcards) else "",
    threads: 8
    conda:
        "../envs/kraken.yaml"
    shell:
        "(kraken2 --db {input.db} --threads {threads} {params.unclassified} "
        " {params.classified} --report {output.report} --gzip-compressed"
        " {params.paired} {input.reads} > {output.kraken_output})"
        "2> {log}"


use rule species_diversity_before_pe as species_diversity_before_se with:
    output:
        kraken_output=temp(
            "results/{date}/species-diversity/se/{sample}/{sample}.kraken"
        ),
        report="results/{date}/species-diversity/se/{sample}/{sample}.uncleaned.kreport2",


# plot Korna charts BEFORE removing human contamination
rule create_krona_chart:
    input:
        kraken_output=get_kraken_output,
        taxonomy_database="resources/krona/",
    output:
        "results/{date}/species-diversity/{sample}/{sample}.html",
    log:
        "logs/{date}/krona/{sample}.log",
    conda:
        "../envs/kraken.yaml"
    shell:
        "ktImportTaxonomy -m 3 -t 5 -tax {input.taxonomy_database} -o {output} {input} 2> {log}"


rule combine_references:
    input:
        "resources/genomes/main.fasta",
        "resources/genomes/human-genome.fna.gz",
    output:
        "resources/genomes/main-and-human-genome.fna.gz",
    log:
        "logs/combine-reference-genomes.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "zcat -f {input} > {output}"


# filter out human contamination
rule extract_reads_of_interest:
    input:
        "results/{date}/mapped/ref~main+human/{sample}.bam",
    output:
        temp("results/{date}/mapped/ref~main+human/nonhuman/{sample}.bam"),
    log:
        "logs/{date}/extract_reads_of_interest/{sample}.log",
    params:
        reference_genome=config["virus-reference-genome"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract-reads-of-interest.py"


rule order_nonhuman_reads:
    input:
        "results/{date}/mapped/ref~main+human/nonhuman/{sample}.bam",
    output:
        fq1=temp("results/{date}/nonhuman-reads/{sample}.1.fastq.gz"),
        fq2=temp("results/{date}/nonhuman-reads/{sample}.2.fastq.gz"),
        bam_sorted=temp("results/{date}/nonhuman-reads/{sample}.sorted.bam"),
    log:
        "logs/{date}/order_nonhuman_reads/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 8
    shell:
        "(samtools sort  -@ {threads} -n {input} -o {output.bam_sorted}; "
        " samtools fastq -@ {threads} {output.bam_sorted} -1 {output.fq1} -2 {output.fq2})"
        " > {log} 2>&1"


# analysis of species diversity present AFTER removing human contamination
rule species_diversity_after:
    input:
        db="resources/minikraken-8GB",
        reads=expand(
            "results/{{date}}/nonhuman-reads/{{sample}}.{read}.fastq.gz", read=[1, 2]
        ),
    output:
        kraken_output=temp(
            "results/{date}/species-diversity-nonhuman/{sample}/{sample}.kraken"
        ),
        report="results/{date}/species-diversity-nonhuman/{sample}/{sample}.cleaned.kreport2",
    log:
        "logs/{date}/kraken/{sample}_nonhuman.log",
    conda:
        "../envs/kraken.yaml"
    threads: 8
    shell:
        "(kraken2 --db {input.db} --threads {threads} --report {output.report} --gzip-compressed "
        "--paired {input.reads} > {output.kraken_output}) 2> {log}"


# plotting Krona charts AFTER removing human contamination
rule create_krona_chart_after:
    input:
        kraken_output="results/{date}/species-diversity-nonhuman/{sample}/{sample}.cleaned.kreport2",
        taxonomy_database="resources/krona/",
    output:
        "results/{date}/species-diversity-nonhuman/{sample}/{sample}.html",
    log:
        "logs/{date}/krona/{sample}_nonhuman.log",
    conda:
        "../envs/kraken.yaml"
    shell:
        "ktImportTaxonomy -m 3 -t 5 -tax {input.taxonomy_database} -o {output} {input} 2> {log}"
