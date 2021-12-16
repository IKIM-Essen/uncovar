rule collect_lineage_candidate_variants:
    input:
        annotation="resources/annotation_known_variants.gff",
        reference="resources/genomes/main.fasta",
    output:
        "resources/lineage-candidate-variants/all.bcf",
    conda:
        "../envs/pysam.yaml"
    log:
        "logs/collect-lineage-candidate-variants.log"
    script:
        "../scripts/collect-lineage-variants.py"


rule annotate_lineage_variants:
    input:
        calls="results/{date}/calls/ref~main/{sample}.lineage-variants.bcf",
        index="results/{date}/calls/ref~main/{sample}.lineage-variants.bcf.csi",
        annotation="resources/lineage-candidate-variants/all_sorted.bcf",
        index2="resources/lineage-candidate-variants/all_sorted.bcf.csi",
    output:
        "results/{date}/lineage-variants/{sample}.bcf",
    log:
        "logs/{date}/annotate-lineage-variants/{sample}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools annotate -a {input.annotation} -c LINEAGES,SIGNATURES {input.calls} > {output} 2> {log}"


# TODO add conda env and log file to this rule
rule generate_lineage_variant_table:
    input:
        lambda wildcards: expand(
            "results/{date}/lineage-variants/{sample}.bcf",
            date=get_dates(),
            sample=get_samples_for_date(get_dates),
            ),
    # input:
    #     "results/{date}/lineage-variants/{sample}.bcf",
    # output:
    #     "results/{data}/lineage-variants/{sample}.tsv",
    # script:
    #     "../scripts/generate-lineage-variant-table.py"


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
