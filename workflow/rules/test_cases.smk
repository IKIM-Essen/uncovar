rule test_non_cov2:
    input:
        pangolin=get_non_cov2_calls(from_caller="pangolin"),
        kallisto=get_non_cov2_calls(from_caller="kallisto"),
        call_plots=expand(
            "results/test-cases/plots/strain-calls/non-cov2-{accession}.strains.{caller}.svg",
            accession=get_non_cov2_accessions(),
            caller=["pangolin", "kallisto"],
        ),
    output:
        "results/test-cases/non-sars-cov-2.csv",
    log:
        "../logs/test-cases/summarize_non_cov2.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize-non-cov2.py"


rule report_non_cov2:
    input:
        "results/test-cases/non-sars-cov-2.csv",
    output:
        report(
            directory("results/test-cases/html"),
            htmlindex="index.html",
            category="Test results",
        ),
    log:
        "../logs/report_non_cov2.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report -s '\t' {input} {output}"
