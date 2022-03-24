# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule vcf_report:
    input:
        ref="resources/genomes/main.fasta",
        index="resources/genomes/main.fasta.fai",
        bams=get_report_input("results/{date}/recal/ref~main/{sample}.bam"),
        bais=get_report_input("results/{date}/recal/ref~main/{sample}.bam.bai"),
        bcfs=get_report_input(
            "results/{date}/filtered-calls/ref~main/{sample}.subclonal.{filter}.{annotation}.bcf"
        ),
    output:
        report(
            directory("results/{date}/vcf-report/{target}.{filter}.{annotation}"),
            htmlindex="index.html",
            caption="../report/variant-calls.rst",
            category="2. Variant Call Details",
            subcategory="Variants with {filter} on {annotation} level",
            labels={
                "sample": "{target}",
                "impact": "{filter}",
                "annotation": "{annotation}",
            },
        ),
    params:
        bcfs=get_report_bcfs,
        bams=get_report_bams,
        format_field="DP AF OBS",
        max_read_depth=config["variant-calling"]["report"]["max-read-depth"],
        js_files="{math} {template}".format(
            math=get_resource("math.min.js"),
            template=get_resource("custom-table-report.js"),
        ),
    log:
        "logs/{date}/vcf-report/{target}.{filter}.{annotation}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-report {input.ref} --bams {params.bams} --vcfs {params.bcfs} --formats {params.format_field} "
        "--infos PROB_* -d {params.max_read_depth} -l {params.js_files} -- {output} 2> {log}"


rule ucsc_vcf:
    input:
        bcfs=get_report_input(
            "results/{date}/filtered-calls/ref~main/{sample}.subclonal.{filter}.{annotation}.bcf"
        ),
        strain_call=(
            "results/{date}/tables/strain-calls/{target}.polished.strains.pangolin.csv"
        ),
    output:
        report(
            "results/{date}/ucsc-vcfs/{target}.{filter}.{annotation}.vcf",
            caption="../report/variant-calls.rst",
            category="5. Variant Call Files",
            subcategory="With {filter} on {annotation}s",
            labels={
                "sample": "{target}",
                "impact": "{filter}",
                "annotation": "{annotation}",
            },
        ),
    log:
        "logs/{date}/ucsc-vcf/{target}.subclonal.{filter}.{annotation}.log",
    conda:
        "../envs/bcftools.yaml"
    script:
        "../scripts/ucsc_vcf.py"


rule aggregate_ucsc_vcfs:
    input:
        expand_samples_for_date(
            "results/{{date}}/ucsc-vcfs/{sample}.{{filter}}.{{annotation}}.vcf"
        ),
    output:
        report(
            "results/{date}/ucsc-vcfs/all.{date}.{filter}.{annotation}.vcf",
            caption="../report/ucsc.rst",
            category="5. Variant Call Files",
            subcategory="Overview",
            labels={"impact": "{filter}", "annotation": "{annotation}"},
        ),
    log:
        "logs/{date}/aggregate_ucsc_vcfs-{filter}-{annotation}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cat {input} > {output} 2> {log}"
