rule uncovar_bcf_2_vcf:
    input:
        "results/{date}/filtered-calls/ref~main/{sample}.subclonal.nofilter.bcf",
    output:
        "results/{date}/filtered-calls/ref~main/{sample}.subclonal.nofilter.vcf",
    log:
        "logs/uncovar_bcf_2_vcf/{date}/{sample}.log",
    conda:
        "../../envs/tools.yaml"
    shell:
        "bcftools view {input} > {output}"


# TODO: check if the other pipelines filter for impact or other criteria
