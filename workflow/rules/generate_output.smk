checkpoint rki_filter:
    input:
        quast_polished_contigs=lambda wildcards: expand(
            "results/{date}/quast-polished/{sample}/report.tsv",
            zip,
            date=[wildcards.date] * len(get_samples_for_date(wildcards.date)),
            sample=get_samples_for_date(wildcards.date),
        ),
        polished_contigs=lambda wildcards: expand(
            "results/{date}/polished-contigs/{sample}.fasta",
            zip,
            date=[wildcards.date] * len(get_samples_for_date(wildcards.date)),
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        temp("results/{date}/rki-filter/{date}.txt"),
    params:
        min_identity=config["RKI-quality-criteria"]["min-identity"],
        max_n=config["RKI-quality-criteria"]["max-n"],
    log:
        "logs/{date}/rki-filter.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rki-filter.py"


rule generate_rki:
    input:
        filtered_samples="results/{date}/rki-filter/{date}.txt",
        polished_contigs=lambda wildcards: expand(
            "results/{date}/polished-contigs/{sample}.fasta",
            zip,
            date=[wildcards.date] * len(get_samples_for_date(wildcards.date)),
            sample=get_samples_for_date(wildcards.date, filtered=True),
        ),
    output:
        fasta="results/rki/{date}_uk-essen_rki.fasta",
        table="results/rki/{date}_uk-essen_rki.csv",
    params:
        min_length=config["rki-output"]["minimum-length"],
    log:
        "logs/{date}/rki-output/{date}.log",
    script:
        "../scripts/generate-rki-output.py"


rule generate_virologist_output:
    input:
        # reads_unfiltered=lambda wildcards: [pep.sample_table.loc[sample][["fq1", "fq2"]] for sample in get_samples_for_date(wildcards.date)],
        reads_unfiltered=lambda wildcards: expand(
            "results/{{date}}/trimmed/{sample}.fastp.json",
            sample=get_samples_for_date(wildcards.date),
        ),
        reads_filtered=lambda wildcards: expand(
            "results/{{date}}/assembly/{sample}/log",
            sample=get_samples_for_date(wildcards.date),
        ),
        initial_contigs=lambda wildcards: expand(
            "results/{{date}}/assembly/{sample}/{sample}.contigs.fa",
            sample=get_samples_for_date(wildcards.date),
        ),
        polished_contigs=lambda wildcards: expand(
            "results/{{date}}/polished-contigs/{sample}.fasta",
            sample=get_samples_for_date(wildcards.date),
        ),
        kraken=lambda wildcards: expand(
            "results/{{date}}/species-diversity/{sample}/{sample}.uncleaned.kreport2",
            sample=get_samples_for_date(wildcards.date),
        ),
        pangolin=lambda wildcards: expand(
            "results/{{date}}/tables/strain-calls/{sample}.strains.pangolin.csv",
            sample=get_samples_for_date(wildcards.date),
        ),
        bcf=lambda wildcards: expand(
            "results/{{date}}/filtered-calls/ref~main/{sample}.subclonal.high+moderate-impact.bcf",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        all_data="results/{date}/virologist/report.csv",
        qc_data="results/{date}/virologist/qc_report.csv",
        var_data="results/{date}/virologist/var_report.csv",
    log:
        "logs/{date}/viro_report.log",
    conda:
        "../envs/pysam.yaml"
    threads: 1
    script:
        "../scripts/generate_virologist_output.py"


rule snakemake_html_report_virologist:
    input:
        qc_data="results/{date}/virologist/qc_report.csv",
        var_data="results/{date}/virologist/var_report.csv",
    output:
        qc_data=report(
            directory("results/{date}/qc_data/"),
            htmlindex="index.html",
            caption="../report/qc-report.rst",
            category="QC report overview",
        ),
        var_data=report(
            directory("results/{date}/var_data/"),
            htmlindex="index.html",
            caption="../report/var-report.rst",
            category="Variant report overview",
        ),
    conda:
        "../envs/rbt.yaml"
    log:
        "logs/{date}/viro_report_html.log",
    shell:
        "(rbt csv-report {input.qc_data} {output.qc_data} && "
        "rbt csv-report {input.var_data} {output.var_data}) > {log} 2>&1"
