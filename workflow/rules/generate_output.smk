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
