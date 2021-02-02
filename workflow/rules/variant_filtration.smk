rule filter_odds:
    input:
        get_filter_odds_input,
    output:
        "results/filtered-calls/ref~{reference}/{sample}.{clonality}.odds.bcf",
    params:
        events=get_target_events,
    log:
        "logs/filter-calls/odds/ref~{reference}/{sample}.{clonality}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls posterior-odds --events {params.events} --odds barely < {input} > {output} 2> {log}"


rule control_fdr:
    input:
        "results/filtered-calls/ref~{reference}/{sample}.{clonality}.odds.bcf",
    output:
        "results/filtered-calls/ref~{reference}/{sample}.{clonality}.{vartype}.fdr-controlled.bcf",
    params:
        fdr=config["variant-calling"]["fdr"],
        events=get_target_events,
    log:
        "logs/control-fdr/ref~{reference}/{sample}.{clonality}.{vartype}.log",
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
        "results/filtered-calls/ref~{reference}/{sample}.{clonality}.bcf",
    log:
        "logs/merge-calls/ref~{reference}/{sample}.{clonality}.log",
    params:
        "-a -Ob",
    wrapper:
        "0.69.0/bio/bcftools/concat"
