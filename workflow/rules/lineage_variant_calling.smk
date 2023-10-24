# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule collect_lineage_candidate_variants:
    input:
        annotation="resources/annotation_known_variants.gff.gz",
        reference="resources/genomes/main.fasta",
    output:
        "resources/lineage-candidate-variants/all.bcf",
    conda:
        "../envs/pysam.yaml"
    log:
        "logs/collect-lineage-candidate-variants.log",
    script:
        "../scripts/collect-lineage-variants.py"


rule annotate_lineage_variants:
    input:
        calls="results/{date}/calls/ref~main/{sample}.lineage-variants.bcf",
        index="results/{date}/calls/ref~main/{sample}.lineage-variants.bcf.csi",
        annotation="resources/lineage-candidate-variants/all.sorted.bcf",
        index2="resources/lineage-candidate-variants/all.sorted.bcf.csi",
    output:
        "results/{date}/lineage-variant-report/{sample}.bcf",
    log:
        "logs/{date}/annotate-lineage-variants/{sample}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools annotate -a {input.annotation} -c LINEAGES,SIGNATURES {input.calls} > {output} 2> {log}"


rule generate_lineage_variant_table:
    input:
        variant_file="results/{date}/lineage-variant-report/{sample}.bcf",
        annotation="resources/annotation_known_variants.gff.gz",
    output:
        variant_table="results/{date}/lineage-variant-report/{sample}.csv",
    log:
        "logs/{date}/variant-table/{sample}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/generate-lineage-variant-table.py"


rule aggregate_lineage_variants:
    input:
        csv=lambda wildcards: expand(
            "results/{{date}}/lineage-variant-report/{sample}.csv",
            sample=get_samples_for_date(wildcards.date),
        ),
        annotation="resources/annotation_known_variants.gff.gz",
    output:
        "results/{date}/lineage-variant-overview/all.csv",
    log:
        "logs/{date}/aggregate-lineage-variants/all.log",
    params:
        sample=lambda wildcards: get_samples_for_date(wildcards.date)
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/aggreagte-lineage-variants.py"


rule get_aggregated_lineage_variant_table:
    input:
        expand(
            "results/{date}/lineage-variant-report/all.csv",
            date="2022-12-21",
        ),


rule render_datavzrd_config:
    input:
        template=workflow.source_path("../../resources/lineage-variant-overview.template.datavzrd.yaml"),
        table="results/{date}/lineage-variant-overview/all.csv",
    output:
        "results/{date}/datavzrd/variant-table-model.yaml",
    log:
        "logs/{date}/yte/render-datavzrd-config/variant-table-model.log",
    template_engine:
        "yte"   


rule render_lineage_variant_table:
    input:
        config="results/{date}/datavzrd/variant-table-model.yaml",
        table="results/{date}/lineage-variant-overview/all.csv",
    output:
        report(
            directory("results/{date}/lineage-variant-report/all"),
            htmlindex="index.html",
            caption="../report/lineage-variant-overview.rst",
            category="2. Variant Call Details",
            subcategory=" VOC variant overview",
        ),
    log:
        "logs/{date}/lineage-variant-overview/all.log",
    wrapper:
        "v2.1.0/utils/datavzrd"



use rule overview_table_html as generate_lineage_variant_report with:
    input:
        "results/{date}/lineage-variant-report/{sample}.csv",
    output:
        report(
            directory(
                "results/{date}/lineage-variant-report/{sample}.lineage-variants"
            ),
            htmlindex="index.html",
            caption="../report/lineage-variant-report.rst",
            category="2. Variant Call Details",
            subcategory=" VOC Similarity",
        ),
    log:
        "logs/{date}/lineage-variant-report/{sample}.log",
    params:
        formatter=get_resource("lineage-variant-table-formatter.js"),
        pin_until="Frequency",
