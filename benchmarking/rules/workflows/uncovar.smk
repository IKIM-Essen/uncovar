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
