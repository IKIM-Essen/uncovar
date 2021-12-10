rule collect_lineage_candidate_variants:
    output:
        "resources/lineage-candidate-variants/all.bcf"
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/collect-lineage-variants.py"


rule generate_lineage_variant_table:
    input:
        "results/{date}/calls/ref~main/{sample}.lineage-variants.bcf"
    output:
        "results/{data}/lineage-variants/{sample}.tsv"
    script:
        "../scripts/generate-lineage-variant-table.py"


use rule overview_table_html as generate_lineage_variant_report with:
    input:
        "results/{data}/lineage-variants/{sample}.tsv"
    output:
        report(
            directory("results/{date}/lineage-variants/{sample}.lineage-variants.tsv"),
            htmlindex="index.html",
            caption="../report/lineage-variant-report.rst",
            category="1. Overview",
            subcategory="1. Report",
        ),
    log:
        "logs/{date}/lineage-variant-report/{sample}.log",