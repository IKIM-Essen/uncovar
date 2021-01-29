rule filter_odds:
    input:
        "results/annotated-calls/{sample}.bcf",
    output:
        "results/filtered-calls/{sample}.odds.bcf",
    log:
        "logs/filter-calls/odds/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls posterior-odds --events PRESENT --odds barely < {input} > {output} 2> {log}"


rule control_fdr:
    input:
        "results/filtered-calls/{sample}.odds.bcf",
    output:
        "results/filtered-calls/{sample}.{vartype}.bcf",
    params:
        fdr=config["variant-calling"]["fdr"],
    log:
        "logs/control-fdr/{sample}.{vartype}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor filter-calls control-fdr {input} --var {wildcards.vartype} "
        "--events PRESENT --fdr {params.fdr} > {output} 2> {log}"


rule bcf_index:
    input:
        "results/filtered-calls/{sample}.{vartype}.bcf",
    output:
        "results/filtered-calls/{sample}.{vartype}.bcf.csi",
    log:
        "logs/bcf-index/{sample}.{vartype}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools index {input} 2> {log}"


rule merge_calls:
    input:
        calls=get_merge_calls_input(".bcf"),
        idx=get_merge_calls_input(".bcf.csi"),
    output:
        "results/filtered-calls/{sample}.bcf",
    log:
        "logs/merge-calls/{sample}.log",
    params:
        "-a -Ob",
    wrapper:
        "0.69.0/bio/bcftools/concat"
