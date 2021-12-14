# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule vembrane_filter:
    input:
        "results/{date}/annotated-calls/ref~main/{sample}.bcf",
    output:
        temp("results/{date}/filtered-calls/ref~main/{sample}.{filter}.bcf"),
    params:
        expression=get_vembrane_expression,
        extra="",
    log:
        "logs/{date}/vembrane/{sample}.{filter}.log",
    wrapper:
        "0.71.1/bio/vembrane/filter"


rule control_fdr:
    input:
        get_control_fdr_input,
    output:
        temp(
            "results/{date}/filtered-calls/ref~{reference}/{sample}.{clonality}.{filter}.{vartype}.fdr-controlled.bcf"
        ),
    params:
        fdr=config["variant-calling"]["fdr"],
        events=get_target_events,
    log:
        "logs/{date}/control-fdr/ref~{reference}/{sample}.{clonality}.{filter}.{vartype}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr --local {input} --var {wildcards.vartype} "
        "--events {params.events} --fdr {params.fdr} > {output} 2> {log}"


rule merge_calls:
    input:
        calls=get_merge_calls_input(".bcf"),
        idx=get_merge_calls_input(".bcf.csi"),
    output:
        "results/{date}/filtered-calls/ref~{reference}/{sample}.{clonality}.{filter}.bcf",
    log:
        "logs/{date}/merge-calls/ref~{reference}/{sample}.{clonality}.{filter}.log",
    params:
        "-a -Ob",
    wrapper:
        "0.69.0/bio/bcftools/concat"
