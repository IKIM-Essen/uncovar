checkpoint extract_strain_genomes_from_gisaid:
    input:
        "resources/gisaid/provision.json",
    output:
        "results/{date}/tables/strain-genomes.txt",
    params:
        save_strains_to=config["strain-calling"]["extracted-strain-genomes"],
    log:
        "logs/{date}/extract-strain-genomes.log",
    params:
        save_strains_to=lambda wildcards: config["strain-calling"][
            "extracted-strain-genomes"
        ],
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/extract-strains-from-gisaid-provision.py"


rule cat_genomes:
    input:
        get_strain_genomes,
    output:
        temp("results/{date}/kallisto/strain-genomes.fasta"),
    log:
        "logs/{date}/cat-genomes.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cat {input} > {output}"


rule kallisto_index:
    input:
        fasta="results/{date}/kallisto/strain-genomes.fasta",
    output:
        index=temp("results/{date}/kallisto/strain-genomes.idx"),
    params:
        extra="",
    log:
        "logs/{date}/kallisto-index.log",
    threads: 8
    wrapper:
        "0.70.0/bio/kallisto/index"


rule kallisto_quant:
    input:
        fastq=get_reads_after_qc,
        index="results/{date}/kallisto/strain-genomes.idx",
    output:
        directory("results/{date}/quant/{sample}"),
    params:
        extra="",
    log:
        "logs/{date}/kallisto_quant/{sample}.log",
    wrapper:
        "0.70.0/bio/kallisto/quant"


rule call_strains_kallisto:
    input:
        quant="results/{date}/quant/{sample}",
        fq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
    output:
        "results/{date}/tables/strain-calls/{sample}.strains.kallisto.tsv",
    log:
        "logs/{date}/call-strains/{sample}.log",
    params:
        min_fraction=config["strain-calling"]["min-fraction"],
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/call-strains.py.ipynb"


rule plot_strains_kallisto:
    input:
        "results/{date}/tables/strain-calls/{sample}.strains.kallisto.tsv",
    output:
        report(
            "results/{date}/plots/strain-calls/{sample}.strains.kallisto.svg",
            caption="../report/strain-calls-kallisto.rst",
            category="1. Overview",
            subcategory="4. Lineage Fraction per Sample",
        ),
    log:
        "logs/{date}/plot-strains-kallisto/{sample}.log",
    params:
        min_fraction=config["strain-calling"]["min-fraction"],
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-strains-kallisto.py.ipynb"


rule plot_all_strains_kallisto:
    input:
        lambda wildcards: expand(
            "results/{{date}}/tables/strain-calls/{sample}.strains.kallisto.tsv",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        report(
            "results/{date}/plots/all.{mode,(major|any)}-strain.strains.kallisto.svg",
            caption="../report/all-strain-calls-kallisto.rst",
            category="1. Overview",
            subcategory="3. Strain Calls",
        ),
    log:
        "logs/{date}/plot-strains/all.{mode}.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-all-strains-kallisto.py.ipynb"


rule pangolin:
    input:
        contigs=get_assemblies_for_submission("single sample"),
        pangoLEARN="results/{date}/pangolin/pangoLEARN",
        lineages="results/{date}/pangolin/lineages",
    output:
        "results/{date}/tables/strain-calls/{sample}.strains.pangolin.csv",
    log:
        "logs/{date}/pangolin/{sample}.log",
    params:
        pango_data_path=lambda w, input: os.path.dirname(input.pangoLEARN),
    conda:
        "../envs/pangolin.yaml"
    threads: 8
    shell:
        "pangolin {input.contigs} --data {params.pango_data_path} --outfile {output} > {log} 2>&1"


rule plot_all_strains_pangolin:
    input:
        lambda wildcards: expand(
            "results/{{date}}/tables/strain-calls/{sample}.strains.pangolin.csv",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        report(
            "results/{date}/plots/all.strains.pangolin.svg",
            caption="../report/all-strain-calls-pangolin.rst",
            category="1. Overview",
            subcategory="3. Strain Calls",
        ),
    log:
        "logs/{date}/plot-strains-pangolin/all.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-all-strains-pangolin.py.ipynb"
