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
include: "workflows/v_pipe.smk"
include: "sanger.smk"
include: "variant_benchmarking.smk"


rule save_workflow_output:
    input:
        get_workflow_output,
        get_output_from_pipline("outdir"),
    output:
        "results/benchmarking/{key}/{tech}/{workflow}/{sample}.fasta",
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


# output:
#     "",
# log:
#     "logs/aggregrate_vc_comparisons.log"
# conda:
#     "../envs/python.yaml"
# script:
#     "../scripts/aggregrate_vc_comparisons.py"
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
# What to compare when benchmarking UnCoVar to other pipelines?
# -> de novo sequences, consensus sequences, variants calls, lineage (pangolin) calls of uncovar vs other pipeline.
# How to compare de novo / consensus sequences?
# -> Alignments
# --> Edit distances to SARS-CoV-2 reference genome by uncovar vs. other pipeline
# --> Identity to SARS-CoV-2 reference genome by uncovar vs. other pipeline
# --> Share N in de novo / consensus sequence by uncovar vs. other pipeline
# --> Visualise difference in seqs. by uncovar vs. other pipeline (b.c. of masking in uncovar)
# How to compare variant calls?
# -> Number of variants, which are also present in sanger sequencing
# -> Number of total found variants by uncovar vs. other pipeline
# -> Number of variants that were found by both pipelines
# -> Uniquely found variants by uncovar vs. other pipeline
# -> Are there variants that are unique to certain vocs, which are only found by one pipeline?
# -> How to integrate probs/ vafs / depth other VCF metrics?
# -> Precision und Recall auf den von Sanger abgedeckten bereichen berechnen
# How to compare Pangolin Calls
# -> Pangolin call by uncovar vs. other pipeline
# -> How to compare Kallisto?
