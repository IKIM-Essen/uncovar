rule rki_filter:
    input:
        quast_polished_contigs=expand(
            "results/quast-polished/{sample}/report.tsv", sample=get_samples()
        ),
        polished_contigs=expand(
            "results/polished-contigs/{sample}.fasta", sample=get_samples()
        ),
    output:
        "results/rki-filter/run.txt",
    params:
        min_identity=config["RKI-quality-criteria"]["min-identity"],
        max_n=config["RKI-quality-criteria"]["max-n"],
    log:
        "logs/rki-filter/run.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rki-filter.py"


rule generate_rki:
    input:
        get_samples_for_date,
    output:
        fasta="results/rki/{date}_uk-essen_rki.fasta",
        table="results/rki/{date}_uk-essen_rki.csv",
    params:
        min_length=config["rki-output"]["minimum-length"],
    log:
        "logs/{date}_rki.log",
    threads: 1
    script:
        "../scripts/generate_rki_output.py"
