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
            subcategory="1. VOC Similarity",
        ),
    log:
        "logs/{date}/lineage-variant-report/{sample}.log",
    params:
        formatter=get_resource("lineage-variant-table-formatter.js"),
        pin_until="Frequency",
