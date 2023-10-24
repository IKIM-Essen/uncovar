rule get_ampliconstats:
    input:
        bed="resources/primer-v4.bed",
        bam="results/{date}/read-sorted/pe~position/{sample}.hardclipped.bam",
    output:
        "results/{date}/ampliconstats/{sample}.ampliconstats.txt",
    log:
        "logs/{date}/samtools/ampliconstats/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools ampliconstats {input.bed} {input.bam} > {output} 2> {log}"


rule plot_amplicon_coverage:
    input:
        bedpe="resources/primer.bedpe",
        amp_stats=lambda wildcards: expand(
            "results/{{date}}/ampliconstats/{sample}.ampliconstats.txt",
            sample=get_samples_for_date(wildcards.date),
        ),
        names="resources/ww-sample-name-dict.csv",
        weeks="resources/ww-sample-week-dict.csv",
    output:
        plot="results/{date}/plots/amplicon-abundance/all.svg",
        stats="results/{date}/tables/amplicon-abundance/all.csv",
    params:
        samples=lambda wildcards: get_samples_for_date(wildcards.date),
    log:
        "logs/{date}/plot-amplicon-coverage.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-amplicon-coverage.py"


rule amplicon_profiles:
    input:
        stats="results/{date}/ampliconstats/{sample}.ampliconstats.txt",
        var_df="results/{date}/lineage-variant-report/{sample}.csv",
    output:
        profiles="results/{date}/amplicon-profiles-sorted/{sample}.txt",
    log:
        "logs/{date}/amplicon-profiles/{sample}.txt",
    threads: 32
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get-amplicon-profiles.py"


rule amplicon_profiles_snv:
    input:
        csv=lambda wildcards: expand(
            "results/{{date}}/lineage-variant-report/{sample}.csv",
            sample=get_samples_for_date(wildcards.date),
            ),
    output:
        ampprofile="results/{date}/amplicon-profiles-snv/all.csv",
    log:
        "logs/{date}/amplicon-profiles-snv/all.log",
    params:
        sample=lambda wildcards: get_samples_for_date(wildcards.date),
    threads: 32
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get-amplicon-profiles-snv.py"


rule get_amplicon_stat_output:
    input:
        expand(
            "results/{date}/amplicon-profiles-snv/all.csv",
            date="2023-09-08",
            sample=get_samples_for_date("2023-09-08"),
            ),
        expand(
            "results/{date}/plots/amplicon-abundance/all.svg",
            date="2023-09-08",
            sample=get_samples_for_date("2023-09-08"),
            ),
        expand(
            "results/{date}/lineage-abundance/all.demix.csv",
            date="2023-09-08",
            sample=get_samples_for_date("2023-09-08"),
            ),