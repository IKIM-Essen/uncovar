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
            caption="../report/strain-calls-kallisto.rst",
            category="Strain calls",
        ),
    log:
        "logs/plot-strains-kallisto/{sample}.log",
    params:
        min_fraction=config["strain-calling"]["min-fraction"],
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-strains-kallisto.py.ipynb"


rule plot_all_strains_kallisto:
    input:
        expand("results/tables/strain-calls/{sample}.strains.tsv", sample=get_samples()),
    output:
        report(
            "results/plots/strain-calls/all.{mode,(major|any)}-strain.strains.kallisto.svg",
            caption="../report/all-strain-calls-kallisto.rst",
            category="Strain calls",
            subcategory="Overview",
        ),
    log:
        "logs/plot-strains/all.{mode}.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-all-strains.py.ipynb"


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


rule plot_strains_pangolin:
    input:
        "results/tables/strain-calls/{sample}.strains.pangolin.csv",
    output:
        report(
            "results/plots/strain-calls/{sample}.strains.pangolin.svg",
            caption="../report/strain-calls-pangolin.rst",
            category="Pangolin strain calls",
        ),
    log:
        "logs/plot-strains-pangolin/{sample}.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-strains-pangolin.py.ipynb"
