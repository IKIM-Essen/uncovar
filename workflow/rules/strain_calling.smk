# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


checkpoint extract_strain_genomes_from_gisaid:
    input:
        "resources/gisaid/provision.json",
    output:
        "results/{date}/tables/strain-genomes.txt",
    params:
        save_strains_to=config["strain-calling"]["extracted-strain-genomes"],
    log:
        "logs/{date}/extract-strain-genomes.log",
    conda:
        "../envs/python.yaml"
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


rule kallisto_metrics:
    input:
        get_reads_after_qc,
    output:
        avg_read_length=temp("results/{date}/tables/avg_read_length/{sample}.txt"),
        standard_deviation=temp("results/{date}/tables/standard_deviation/{sample}.txt"),
    log:
        "logs/{date}/kallisto/metrics/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "awk 'BEGIN {{ t=0.0;sq=0.0; n=0; }} ;NR%4==2 {{ n++;L=length($0);t+=L;sq+=L*L; }}END{{ m=t/n;printf(\"%f\\n\",m) ; }}' {input} && "
        "awk 'BEGIN {{ t=0.0;sq=0.0; n=0; }} ;NR%4==2 {{ n++;L=length($0);t+=L;sq+=L*L; }}END{{ m=t/n;printf(\"%f\\n\",sq/n-m*m) ; }}' {input} && "
        "awk 'BEGIN {{ t=0.0;sq=0.0; n=0; }} ;NR%4==2 {{ n++;L=length($0);t+=L;sq+=L*L; }}END{{ m=t/n;printf(\"%f\\n\",m) ; }}' {input} > {output.avg_read_length} && "
        "awk 'BEGIN {{ t=0.0;sq=0.0; n=0; }} ;NR%4==2 {{ n++;L=length($0);t+=L;sq+=L*L; }}END{{ m=t/n;printf(\"%f\\n\",sq/n-m*m) ; }}' {input} > {output.standard_deviation}"


rule kallisto_quant:
    input:
        unpack(get_kallisto_quant_input),
    output:
        directory("results/{date}/quant/{sample}"),
    params:
        extra=lambda w, input: get_kallisto_quant_extra(w, input),
    log:
        "logs/{date}/kallisto_quant/{sample}.log",
    threads: 8
    wrapper:
        "0.70.0/bio/kallisto/quant"


rule kallisto_call_strains:
    input:
        quant="results/{date}/quant/{sample}",
        fq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
    output:
        "results/{date}/tables/strain-calls/{sample}.strains.kallisto.tsv",
    log:
        "logs/{date}/call-strains/{sample}.log",
    params:
        min_fraction=config["strain-calling"]["min-fraction"],
    resources:
        notebooks=1,
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/call-strains.py.ipynb"


rule kallisto_plot_strains:
    input:
        "results/{date}/tables/strain-calls/{sample}.strains.kallisto.tsv",
    output:
        report(
            "results/{date}/plots/strain-calls/{sample}.strains.kallisto.svg",
            caption="../report/strain-calls-kallisto.rst",
            category="1. Overview",
            subcategory="4. Lineage Fraction per Sample",
            labels={"sample": "{sample}"},
        ),
    log:
        "logs/{date}/plot-strains-kallisto/{sample}.log",
    params:
        min_fraction=config["strain-calling"]["min-fraction"],
    resources:
        notebooks=1,
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-strains-kallisto.py.ipynb"


rule kallisto_plot_all_strains:
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
            subcategory="3. Lineage Calls",
        ),
    log:
        "logs/{date}/plot-strains/all.{mode}.log",
    resources:
        notebooks=1,
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-all-strains-kallisto.py.ipynb"


rule pangolin_call_strains:
    input:
        contigs=get_pangolin_input,
        pangoLEARN="results/{date}/pangolin/pangoLEARN",
        lineages="results/{date}/pangolin/lineages",
    output:
        "results/{date}/tables/strain-calls/{sample}.{stage}.strains.pangolin.csv",
    log:
        "logs/{date}/pangolin/{sample}.{stage}.log",
    params:
        pango_data_path=lambda w, input: os.path.dirname(input.pangoLEARN),
    conda:
        "../envs/pangolin.yaml"
    threads: 8
    shell:
        "pangolin {input.contigs} --data {params.pango_data_path} --outfile {output} > {log} 2>&1"


rule pangolin_plot_all_strains:
    input:
        lambda wildcards: expand(
            "results/{{date}}/tables/strain-calls/{sample}.polished.strains.pangolin.csv",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        report(
            "results/{date}/plots/all.strains.pangolin.svg",
            caption="../report/all-strain-calls-pangolin.rst",
            category="1. Overview",
            subcategory="3. Lineage Calls",
        ),
    log:
        "logs/{date}/plot-strains-pangolin/all.log",
    resources:
        notebooks=1,
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-all-strains-pangolin.py.ipynb"
