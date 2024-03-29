# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule vembrane_filter:
    input:
        vcf="results/{date}/annotated-calls/ref~main/annot~{annotation}/{sample}.bcf",
    output:
        vcf=temp(
            "results/{date}/filtered-calls/ref~main/annot~{annotation}/{sample}.{filter}.bcf"
        ),
    params:
        expression=get_vembrane_expression,
        extra="",
    log:
        "logs/{date}/vembrane/{sample}.{filter}.{annotation}.log",
    wrapper:
        "v1.15.1/bio/vembrane/filter"


rule control_fdr:
    input:
        get_control_fdr_input,
    output:
        temp(
            "results/{date}/filtered-calls/ref~{reference}/{sample}.{clonality}.{filter}.{vartype}.{annotation}.fdr-controlled.bcf"
        ),
    params:
        fdr=config["variant-calling"]["fdr"],
        events=get_target_events,
    log:
        "logs/{date}/control-fdr/ref~{reference}/{sample}.{clonality}.{filter}.{vartype}.{annotation}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr --mode local-smart {input} --var {wildcards.vartype} "
        "--events {params.events} --fdr {params.fdr} > {output} 2> {log}"


rule merge_calls:
    input:
        calls=get_merge_calls_input(".bcf"),
        idx=get_merge_calls_input(".bcf.csi"),
    output:
        "results/{date}/filtered-calls/ref~{reference}/{sample}.{clonality}.{filter}.{annotation}.bcf",
    log:
        "logs/{date}/merge-calls/ref~{reference}/{sample}.{clonality}.{filter}.{annotation}.log",
    params:
        extra="-a",
    wrapper:
        "v1.15.1/bio/bcftools/concat"
