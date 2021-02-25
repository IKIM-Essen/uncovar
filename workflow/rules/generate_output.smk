checkpoint rki_filter:
    input:
        quast_polished_contigs=lambda wildcards: expand(
            "results/quast-polished/{sample}/report.tsv",
            sample=get_samples_for_date(wildcards.date),
        ),
        polished_contigs=lambda wildcards: expand(
            "results/polished-contigs/{sample}.fasta",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        "results/rki-filter/{date}.txt",
    params:
        min_identity=config["RKI-quality-criteria"]["min-identity"],
        max_n=config["RKI-quality-criteria"]["max-n"],
    log:
        "logs/rki-filter/{date}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rki-filter.py"


rule generate_rki:
    input:
        lambda wildcards: expand(
            "results/polished-contigs/{sample}.fasta",
            sample=get_samples_for_date(wildcards.date, filtered=True),
        ),
    output:
        fasta="results/rki/{date}_uk-essen_rki.fasta",
        table="results/rki/{date}_uk-essen_rki.csv",
    params:
        min_length=config["rki-output"]["minimum-length"],
    log:
        "logs/{date}_rki.log",
    script:
        "../scripts/generate_rki_output.py"
