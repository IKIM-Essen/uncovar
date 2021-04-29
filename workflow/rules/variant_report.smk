rule vcf_report:
    input:
        ref="resources/genomes/main.fasta",
        bams=get_report_input("results/{date}/recal/ref~main/{sample}.bam"),
        bais=get_report_input("results/{date}/recal/ref~main/{sample}.bam.bai"),
        bcfs=get_report_input(
            "results/{date}/filtered-calls/ref~main/{sample}.subclonal.{filter}.bcf"
        ),
    output:
        report(
            directory("results/{date}/vcf-report/{target}.{filter}"),
            htmlindex="index.html",
            caption="../report/variant-calls.rst",
            category="2. Variant Call Details",
            subcategory="{filter}",
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
        "logs/{date}/vcf-report/{target}.{filter}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-report {input.ref} --bams {params.bams} --vcfs {params.bcfs} --formats {params.format_field} "
        "--infos PROB_* -d {params.max_read_depth} -l {params.js_files} -- {output} 2> {log}"


rule ucsc_vcf:
    input:
        bcfs=get_report_input(
            "results/{date}/filtered-calls/ref~main/{sample}.subclonal.{filter}.bcf"
        ),
        strain_call="results/{date}/tables/strain-calls/{target}.strains.pangolin.csv",
    output:
        report(
            "results/{date}/ucsc-vcfs/{target}.{filter}.vcf",
            caption="../report/variant-calls.rst",
            category="6. Variant Call Files",
            subcategory="{filter}",
        ),
    log:
        "logs/{date}/ucsc-vcf/{target}.subclonal.{filter}.log",
    conda:
        "../envs/bcftools.yaml"
    script:
        "../scripts/ucsc_vcf.py"


rule aggregate_ucsc_vcfs:
    input:
        lambda wildcards: expand(
            "results/{{date}}/ucsc-vcfs/{sample}.{{filter}}.vcf",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        report(
            "results/{date}/ucsc-vcfs/all.{date}.{filter}.vcf",
            caption="../report/ucsc.rst",
            category="6. Variant Call Files",
            subcategory="Overview",
        ),
    log:
        "logs/{date}/aggregate_ucsc_vcfs-{filter}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cat {input} > {output} 2> {log}"
