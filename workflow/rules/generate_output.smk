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
        # initial_contigs=lambda wildcards: expand(
        #     "results/{{date}}/assembly/{sample}/{sample}.contigs.fa",
        #     sample=get_samples_for_date(wildcards.date),
        # ),
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
        "results/{date}/virologist/report.csv",
    log:
        "logs/{date}/viro_report.log",
    conda:
        "../envs/pysam.yaml"
    threads: 1
    script:
        "../scripts/generate_virologist_output.py"


rule report_virologist:
    input:
        "results/{date}/virologist/report.csv",
    output:
        report(
            directory("results/{date}/virologist-report"),
            htmlindex="index.html",
            #caption="../report/virologist-report.rst",
            category="Virologist report",
        ),
    conda:
        "../envs/rbt.yaml"
    log:
        "logs/{date}/viro_report_html.log",
    shell:
        "(rbt csv-report -s ',' {input} {output}) > {log} 2>&1"