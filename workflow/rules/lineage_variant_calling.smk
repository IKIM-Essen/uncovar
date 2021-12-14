rule collect_lineage_candidate_variants:
    input:
        annotation="resources/annotation.gff",
        reference="resources/genomes/main.fasta",
    output:
        "resources/lineage-candidate-variants/all.bcf",
    conda:
        "../envs/pysam.yaml"
    log:
        "logs/"
    script:
        "../scripts/collect-lineage-variants.py"


rule annotate_lineage_variants:
    input:
        calls="results/{date}/calls/ref~main/{sample}.lineage-variants.bcf",
        annotation="resources/lineage-candidate-variants/all.bcf",
    output:
        "results/{date}/lineage-variants/{sample}.bcf",
    log:
        "logs/annotate-lineage-variants/{sample}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools annotate -a {input.annotation} -c LINEAGES,SIGNATURES {input.calls} > {output} 2> {log}"


# TODO add conda env and log file to this rule
rule generate_lineage_variant_table:
    input:
        "results/{date}/lineage-variants/{sample}.bcf",
    output:
        "results/{data}/lineage-variants/{sample}.tsv",
    script:
        "../scripts/generate-lineage-variant-table.py"


use rule overview_table_html as generate_lineage_variant_report with:
    input:
        "results/{data}/lineage-variants/{sample}.tsv",
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
