include: "common.smk"
include: "ref.smk"
include: "artic.smk"


# include: "ncov2019_artic_nf.smk"
# include: "nf_core_viralrecon.smk"
# include: "porecov.smk"
# include: "v_pipe.smk"
# include: "covpipe.smk"


rule agg_vcf:
    input:
        lambda w: expand(
            [
                "results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.merged.vcf",
                "results/benchmarking/artic/minion/medaka/{sample}/{sample}.merged.vcf",
            ],
            sample=get_nanopore_samples(w),
        ),


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
#     threads: 8
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
