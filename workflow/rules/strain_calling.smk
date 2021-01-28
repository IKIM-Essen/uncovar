rule cat_genomes:
    input:
        get_strain_genomes,
    output:
        "resources/strain-genomes.fasta",
    log:
        "logs/cat-genomes.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cat {input} > {output}"


rule kallisto_index:
    input:
        fasta="resources/strain-genomes.fasta",
    output:
        index="resources/strain-genomes.idx",
    params:
        extra="",
    log:
        "logs/kallisto-index.log",
    threads: 8
    wrapper:
        "0.70.0/bio/kallisto/index"


rule kallisto_quant:
    input:
        fastq=expand("results/trimmed/{{sample}}.{read}.fastq.gz", read=[1, 2]),
        index="resources/strain-genomes.idx",
    output:
        directory("results/quant/{sample}"),
    params:
        extra="",
    log:
        "logs/kallisto_quant/{sample}.log",
    threads: 1
    wrapper:
        "0.70.0/bio/kallisto/quant"


rule call_strains_kallisto:
    input:
        "results/quant/{sample}",
    output:
        "results/tables/strain-calls/{sample}.strains.kallisto.tsv",
    log:
        "logs/call-strains/{sample}.log",
    params:
        min_fraction=config["strain-calling"]["min-fraction"],
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/call-strains.py.ipynb"


rule plot_strains_kallisto:
    input:
        "results/tables/strain-calls/{sample}.strains.kallisto.tsv",
    output:
        report(
            "results/plots/strain-calls/{sample}.strains.kallisto.svg",
            caption="../report/strain-calls.rst",
            category="Strain calls",
            subcategory="Per sample",
        ),
    log:
        "logs/plot-strains/{sample}.log",
    params:
        min_fraction=config["strain-calling"]["min-fraction"],
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-strains.py.ipynb"



rule pangolin:
    input:
        "results/assembly/{sample}/final.contigs.fa",
    output:
        "results/tables/strain-calls/{sample}.strains.pangolin.csv",
    log:
        "logs/pangolin/{sample}.log",
    threads: 8
    conda:
        "../envs/pangolin.yaml"
    shell:
        "pangolin {input} --outfile {output}"

rule plot_all_strains:
    input:
        expand("results/tables/strain-calls/{sample}.strains.tsv", sample=get_samples()),
    output:
        report(
            "results/plots/strain-calls/all.{mode,(major|any)}-strain.strains.svg",
            caption="../report/all-strain-calls.rst",
            category="Strain calls",
            subcategory="Overview",
        ),
    log:
        "logs/plot-strains/all.{mode}.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-all-strains.py.ipynb"


# TODO
# 1. the entrez rule get_genome also downloads partial reads of covid genomes (e.g. MW368461) or empty genome files (e.g. MW454604). Add rule to exiculde those rules "faulty" covid genomes

