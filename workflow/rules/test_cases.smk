rule test_non_cov2:
    input:
        pangolin = get_non_cov2_calls(from_caller="pangolin"),
        kallisto = get_non_cov2_calls(from_caller="kallisto"),
        pangolin_plots= expand("results/test-cases/plots/strain-calls/non-cov2-{accession}.strains.pangolin.svg", accession=get_non_cov2_accessions()),
        kallisto_plots= expand("results/test-cases/plots/strain-calls/non-cov2-{accession}.strains.kallisto.svg", accession=get_non_cov2_accessions())
    output:
        report(
            "results/test-cases/non-sars-cov-2.csv",
            category="Test results of non-SARS-CoV-2 genomes",
        )
    log:
        "../logs/test-cases/summarize_non_cov2.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize-non-cov2.py"