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


rule snakemake_reports:
    input:
        lambda wildcards: expand(
            "results/{{date}}/polished-contigs/{sample}.fasta",
            sample=get_samples_for_date(wildcards.date),
        ),
        lambda wildcards: expand(
            "results/{{date}}/plots/strain-calls/{sample}.strains.kallisto.svg",
            sample=get_samples_for_date(wildcards.date),
        ),
        expand(
            "results/{{date}}/plots/all.{mode}-strain.strains.kallisto.svg",
            mode=["major", "any"],
        ),
        lambda wildcards: expand(
            "results/{{date}}/plots/strain-calls/{sample}.strains.pangolin.svg",
            sample=get_samples_for_date(wildcards.date),
        ),
        "results/{date}/plots/all.strains.pangolin.svg",
        lambda wildcards: expand(
            "results/{{date}}/scenarios/{sample}.yaml",
            sample=get_samples_for_date(wildcards.date),
        ),
        lambda wildcards: expand(
            "results/{{date}}/vcf-report/{target}.{filter}",
            target=get_samples_for_date(wildcards.date) + ["all"],
            filter=config["variant-calling"]["filters"],
        ),
    output:
        "results/reports/{date}.zip",
    params:
        for_testing = "--snakefile ../workflow/Snakefile" if config.get("benchmark-genomes", []) else ""
    log:
        "../logs/snakemake_reports/{date}.log",
    shell:
        "snakemake {input} --report {output} {params.for_testing}"
