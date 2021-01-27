rule vcf_report:
    input:
        ref="resources/genomes/main.fasta",
        bams=get_report_input("results/recal/{sample}.bam"),
        bais=get_report_input("results/recal/{sample}.bam.bai"),
        bcfs=get_report_input("results/filtered-calls/{sample}.bcf"),
    output:
        report(
            directory("results/vcf-report/{target}"),
            htmlindex="index.html",
            caption="../report/variant-calls.rst",
            category="Variant calls",
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
        "logs/vcf-report/{target}.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-report {input.ref} --bams {params.bams} --vcfs {params.bcfs} --formats {params.format_field} "
        "--infos PROB_* -d {params.max_read_depth} -l {params.js_files} -- {output} 2> {log}"
