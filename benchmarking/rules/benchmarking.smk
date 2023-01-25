include: "common.smk"
include: "ref.smk"
include: "workflows/artic.smk"
include: "workflows/covpipe.smk"
include: "workflows/havoc.smk"
include: "workflows/ncov2019_artic_nf.smk"
include: "workflows/nf_core_viralrecon.smk"
include: "workflows/porecov.smk"
include: "workflows/signal.smk"
include: "workflows/snakelines.smk"
include: "workflows/uncovar.smk"
include: "workflows/v_pipe.smk"
include: "sanger.smk"
include: "benchmarking_variants.smk"
include: "benchmarking_sequences.smk"
include: "benchmarking_lineages.smk"
include: "benchmarking_time.smk"


rule save_workflow_output:
    input:
        get_workflow_output,
        get_output_from_pipline("outdir"),
    output:
        "results/benchmarking/backups/{key}/{tech}/{workflow}/{sample}.some.extension",
    log:
        "logs/save/{key}/{tech}/{workflow}/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cp {input[0]} {output} 2> {log}"


rule save_all_workflow_outputs:
    input:
        get_all_outputs,
    output:
        touch("results/benchmarking/.saved"),


###################
# Tried workflows #
###################
# source: https://github.com/niemasd/ViReflow
# ViReflow is not runable
# -> needs s3 bucket
# rule ViReflow:
#     input:
#         fq=get_fastqs,
#         script="resources/benchmarking/ViReflow/ViReflow.py",
#         reference="resources/genomes/main.fasta",
#         gff="resources/annotation.gff.gz",
#         bed="resources/primers.bed",
#     output:
#         directory("results/benchmarking/ViReflow/{sample}"),
#     log:
#         "logs/ViReflow/{sample}.log",
#     threads: 4
#     conda:
#         "../../envs/python.yaml"
#     shell:
#         "./{input.script} --destination {output} --reference_fasta {input.reference} "
#         "--reference_gff {input.gff} --primer_bed {input.bed} --output {output} "
#         "--threads {threads} --optional_pangolin true  {input.fq[0]}"
###################
# C-View is not installable or runable.
# -> Needs sudo
# -> paths to softwaredir and anaconda dir sometimes hardcoded
# -> Was not able to start
# rule C_VIEW_install:
#     input:
#         "resources/benchmarking/C-VIEW/install.sh",
#     output:
#         directory("resources/benchmarking/C-VIEW/softwaredir"),
#     log:
#         "logs/C_VIEW_install.log"
#     conda:
#         "../../envs/c-view.yaml"
#     shell:
#         "./{input} {output} $(which anaconda | sed 's/\/bin.*//g')"
