rule uncovar_bcf_2_vcf:
    input:
        "results/{date}/filtered-calls/ref~main/{sample}.subclonal.high+moderate-impact.bcf",
    output:
        "results/{date}/filtered-calls/ref~main/{sample}.subclonal.high+moderate-impact.vcf",
    log:
        "logs/uncovar_bcf_2_vcf/{date}/{sample}.log",
    conda:
        "../../envs/tools.yaml"
    shell:
        "bcftools view {input} > {output}"


# TODO: check if the other pipelines filter for impact or other criteria


rule uncovar_control_fdr:
    input:
        "results/{date}/annotated-calls/ref~main/{sample}.bcf",
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
