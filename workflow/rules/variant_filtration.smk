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


rule filter_odds:
    input:
        get_filter_odds_input,
    output:
        temp(
            "results/{date}/filtered-calls/ref~{reference}/{sample}.{clonality}.{filter}.odds.bcf"
        ),
    params:
        events=get_target_events,
    log:
        "logs/{date}/filter-calls/odds/ref~{reference}/{sample}.{clonality}.{filter}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls posterior-odds --events {params.events} --odds barely < {input} > {output} 2> {log}"


rule control_fdr:
    input:
        "results/{date}/filtered-calls/ref~{reference}/{sample}.{clonality}.{filter}.odds.bcf",
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
        "varlociraptor filter-calls control-fdr {input} --var {wildcards.vartype} "
        "--events {params.events} --fdr {params.fdr} > {output} 2> {log}"


rule merge_calls:
    input:
        calls=get_merge_calls_input(".bcf"),
        idx=get_merge_calls_input(".bcf.csi"),
    output:
        temp(
            "results/{date}/filtered-calls/ref~{reference}/{sample}.{clonality}.{filter}.bcf"
        ),
    log:
        "logs/{date}/merge-calls/ref~{reference}/{sample}.{clonality}.{filter}.log",
    params:
        "-a -Ob",
    wrapper:
        "0.69.0/bio/bcftools/concat"
