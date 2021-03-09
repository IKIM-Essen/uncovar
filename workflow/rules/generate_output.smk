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
        fasta=report(
            "results/rki/{date}_uk-essen_rki.fasta",
            category="4. RKI Submission",
            caption="../report/rki-submission-fasta.rst",
        ),
        table=report(
            "results/rki/{date}_uk-essen_rki.csv",
            category="4. RKI Submission",
            caption="../report/rki-submission-csv.rst",
        ),
    params:
        min_length=config["rki-output"]["minimum-length"],
    log:
        "logs/{date}/rki-output/{date}.log",
    script:
        "../scripts/generate-rki-output.py"


rule snakemake_report:
    input:
        "results/{date}/plots/all.strains.pangolin.svg",
        lambda wildcards: expand(
            "results/{{date}}/vcf-report/{target}.{filter}",
            target=get_samples_for_date(wildcards.date) + ["all"],
            filter=config["variant-calling"]["filters"],
        ),
        "results/{date}/qc/laboratory/multiqc.html",
        "results/rki/{date}_uk-essen_rki.csv",
        "results/rki/{date}_uk-essen_rki.fasta",
        lambda wildcards: expand(
            "results/{{date}}/polished-contigs/{sample}.fasta",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        "results/reports/{date}.zip",
    params:
        for_testing=(
            "--snakefile ../workflow/Snakefile --nolock"
            if config.get("benchmark-genomes", [])
            else ""
        ),
    log:
        "../logs/snakemake_reports/{date}.log",
    shell:
        "snakemake --report-stylesheet resources/custom-stylesheet.css {input} --report {output} {params.for_testing}"
