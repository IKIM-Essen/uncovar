include: "common.smk"
include: "ref.smk"
include: "artic.smk"
include: "covpipe.smk"
include: "ncov2019_artic_nf.smk"
include: "nf_core_viralrecon.smk"
include: "porecov.smk"
# include: "signal.smk"
include: "v_pipe.smk"
include: "sanger.smk"


# rule extract_vcf:
#     input:
#         "results/benchmarking/{infix}.vcf.gz",
#     output:
#         "results/benchmarking/{infix}.vcf",
#     log:
#         "logs/extract_vcf/{infix}.log",
#     conda:
#         "../envs/unix.yaml"
#     shell:
#         "bgzip k {input}"


# rule compress_vcf:
#     input:
#         "results/benchmarking/{infix}.vcf",
#     output:
#         "results/benchmarking/{infix}.vcf.gz",
#     log:
#         "logs/compress_vcf/{infix}.log"
#     conda:
#         "../envs/unix.yaml"
#     shell:
#         "gzip -dk {input}"


rule agg_vcf:
    input:
        nanopore_artic_medaka=lambda w: expand(
            "results/benchmarking/artic/minion/medaka/{sample}/{sample}.merged.vcf",
            sample=get_nanopore_samples(w),
        ),
        nanopore_artic_nanopolish=lambda w: expand(
            "results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.merged.vcf",
            sample=get_nanopore_samples(w),
        ),
        nanopore_ncov2019_artic_nf_medaka=lambda w: expand(
            "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{sample}-{barcode}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/{sample}_{barcode}.merged.vcf",
            zip,
            sample=get_nanopore_samples(w),
            barcode=get_barcodes(w),
        ),
        nanopore_ncov2019_artic_nf_nanopolish=lambda w: expand(
            "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{sample}-{barcode}/articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish/{sample}_{barcode}.merged.vcf",
            zip,
            sample=get_nanopore_samples(w),
            barcode=get_barcodes(w),
        ),
        nanopore_nf_core_viralrecon_nanopolish=lambda w: expand(
            "results/benchmarking/nf-core-viralrecon/nanopore/nanopolish/{sample}/nanopolish/{sample}.merged.vcf",
            sample=get_nanopore_samples(w),
        ),
        nanopore_nf_core_viralrecon_medaka=lambda w: expand(
            "results/benchmarking/nf-core-viralrecon/nanopore/medaka/{sample}/medaka/{sample}.merged.vcf",
            sample=get_nanopore_samples(w),
        ),
        illumina_covpipe=lambda w: expand(
            "results/benchmarking/CovPipe/{sample}-{covpipe_name}/results/intermediate_data/04_variant_calling/{covpipe_name}/{covpipe_name}.vcf",
            zip,
            sample=get_illumina_samples(w),
            covpipe_name=get_covpipe_names(w),
        ),
        illumina_ncov2019_artic_nf=lambda w: expand(
            "results/benchmarking/ncov2019_artic_nf/illumina/{sample}/ncovIllumina_sequenceAnalysis_callVariants/{sample}.variants.tsv",
            sample=get_illumina_samples(w),
        ),
        illumina_nf_core_viralrecon=lambda w: expand(
            "results/benchmarking/nf-core-viralrecon/illumina/{sample}/variants/bcftools/{sample}.vcf",
            sample=get_illumina_samples(w),
        ),
        illumina_v_pipe=lambda w: expand(
            "results/benchmarking/v-pipe/{sample}/work/samples/{sample}/20200102/variants/SNVs/snvs.vcf",
            sample=get_illumina_samples(w),
        ),


rule compare_vc_of_pipelines:
    input:
        fist_pipeline=lambda w: get_vcf_of_pipeline(w.pipeline_one, w),
        second_pipeline=lambda w: get_vcf_of_pipeline(w.pipeline_two, w),
        other="results/2022-01-19/filtered-calls/ref~main/28998_illumina.subclonal.high+moderate-impact.bcf.csi",
    output:
        "results/benchmarking/tables/vc-comparison-{sample}~{pipeline_one}-vs-{pipeline_two}.tsv",
    log:
        "logs/compare_vc_of_pipelines/sample-{sample}~{pipeline_one}-vs-{pipeline_two}.log",
    params:
        fist_pipeline=lambda w: f"{w.pipeline_one}",
        second_pipeline=lambda w: f"{w.pipeline_two}",
    conda:
        "../envs/vcftools.yaml"
    shell:
        "bcftools stats {input.fist_pipeline} {input.second_pipeline} > {output}"


rule aggregrate_vc_comparisons:
    input:
        "results/benchmarking/tables/vc-comparison-28998_illumina~uncovar-vs-nf-core-viralrecon.tsv",


# rule preprocess_variants:
#     input:
#         ##vcf/bcf
#         variants=lambda w: get_vcf_of_pipeline(w.pipeline, w)
#     output:
#         "results/benchmarking/happy/normalized/{pipeline}/{sample}.bcf"
#     params:
#         ## path to reference genome
#         genome="resources/genomes/main.fasta",
#         ## parameters such as -L to left-align variants
#         extra="-L"
#     log: "logs/happy/preprocess/{pipeline}/{sample}.log"
#     threads: 2
#     wrapper:
#         "v1.0.0/bio/hap.py/pre.py"


rule benchmark_variants:
    input:
        truth=lambda w: get_vcf_of_pipeline(w.pipeline_1, w),  # sanger vcf
        query=lambda w: get_vcf_of_pipeline(w.pipeline_2, w),  # variant calls
        truth_regions="results/benchmarking/sanger/aligned/{sample}.bed",  # sanger bed
        genome="resources/genomes/main.fasta",
        genome_index="resources/genomes/main.fasta.fai",
    output:
        multiext(
            "results/benchmarking/happy/happy/truth~{pipeline_1}-vs-{pipeline_2}/{sample}/{sample}",
            ".runinfo.json",
            ".vcf.gz",
            ".summary.csv",
            ".extended.csv",
            ".metrics.json.gz",
            ".roc.all.csv.gz",
            ".roc.Locations.INDEL.csv.gz",
            ".roc.Locations.INDEL.PASS.csv.gz",
            ".roc.Locations.SNP.csv.gz",
            ".roc.tsv",
        ),
    params:
        engine="vcfeval",
        prefix=lambda wc, input, output: output[0].split(".")[0],
        ## parameters such as -L to left-align variants
        extra="--verbose",
    log:
        "logs/happy/truth~{pipeline_1}-vs-{pipeline_2}/{sample}.log",
    threads: 2
    wrapper:
        "v1.0.0/bio/hap.py/hap.py"


rule agg_happy:
    input:
        lambda w: expand(
            "results/benchmarking/happy/happy/truth~sanger-vs-{pipeline}/{sample}/{sample}.runinfo.json",
            pipeline=PIPELINES["nanopore"],
            sample=get_nanopore_samples(w),
        ),
        lambda w: expand(
            "results/benchmarking/happy/happy/truth~sanger-vs-{pipeline}/{sample}/{sample}.runinfo.json",
            pipeline=PIPELINES["illumina"],
            sample=get_illumina_samples(w),
        ),


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
