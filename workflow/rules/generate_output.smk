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
