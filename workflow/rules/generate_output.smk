rule generate_rki:
    input:
        expand(
            "results/polished-contigs/{sample}.fasta", sample=get_samples_latest_run(),
        ),
    output:
        fasta=expand("results/rki/{date}_uk-essen_rki.fasta", date=get_latest_run_date()),
        table=expand("results/rki/{date}_uk-essen_rki.csv", date=get_latest_run_date()),
    log:
        "logs/rki.log"
    threads: 1
    script:
        "/home/alex/repo/snakemake-workflow-sars-cov2/processing_scripts/generate_rki_output.py"
