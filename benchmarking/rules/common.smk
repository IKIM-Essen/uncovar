# PIPELINES = {"nanopore": ["uncovar"], "illumina": ["signal"]}
PIPELINES = {
    "nanopore": [
        "artic-medaka",
        "artic-nanopolish",
        "ncov2019-artic-nf-medaka",
        "ncov2019-artic-nf-nanopolish",
        "nf-core-viralrecon-nanopolish",
        "nf-core-viralrecon-medaka",
        "uncovar",
    ],
    "illumina": [
        "ncov2019-artic-nf",
        "nf-core-viralrecon",
        "v-pipe",
        "havoc",
        "covpipe",
        "signal",
        "snakelines",
        "uncovar",
    ],
}


def get_fastq_pass_path_barcode(wildcards, sample=None):
    if sample is not None:
        return pep.sample_table.loc[sample]["fastq_pass"]
    return pep.sample_table.loc[wildcards.sample]["fastq_pass"]


def get_fast5_pass_path_barcode(wildcards):
    return pep.sample_table.loc[wildcards.sample]["fast5_pass"]


def get_seq_summary(wildcards):
    return pep.sample_table.loc[wildcards.sample]["seq_summary"]


def get_barcode_for_viralrecon_nanopore_sample(wildcards):
    barcode = os.path.basename(os.path.normpath(get_fastq_pass_path_barcode(wildcards)))
    barcode = barcode.replace("barcode0", "")
    barcode = barcode.replace("barcode", "")
    return f"sample,barcode\n{wildcards.sample},{barcode}"


def get_barcode_for_viralrecon_illumina_sample(wildcards):
    fq1, fq2 = get_fastqs(wildcards)
    return f"sample,fastq_1,fastq_2\n{wildcards.sample},{fq1},{fq2}"


def get_fastq_or_fast5(wildcards):
    if wildcards.folder == "fastq_pass":
        return get_fastq_pass_path_barcode(wildcards)
    if wildcards.folder == "fast5_pass":
        return get_fast5_pass_path_barcode(wildcards)


def get_barcode(wildcards):
    return os.path.basename(os.path.normpath(get_fastq_pass_path_barcode(wildcards)))


def get_barcodes(wildcards):
    fastq_paths = [
        get_fastq_pass_path_barcode(None, sample)
        for sample in get_nanopore_samples(wildcards)
    ]
    return [os.path.basename(os.path.normpath(sample)) for sample in fastq_paths]


def get_covpipe_names(wildcards):
    return [sample.replace("_", "__") for sample in get_illumina_samples(wildcards)]


def get_covpipe_name_for_sample(wildcards):
    return wildcards.sample.replace("_", "__")


def get_nanopore_samples(wildcards):
    return pep.sample_table.loc[
        pep.sample_table["technology"] == "ont", "sample_name"
    ].values


def get_illumina_samples(wildcards):
    return pep.sample_table.loc[
        pep.sample_table["technology"] == "illumina", "sample_name"
    ].values


def get_date_for_sample(wildcards):
    return pep.sample_table.loc[wildcards.sample]["date"]


def get_vcf_of_workflow(pipeline, wildcards):
    if pipeline == "artic-medaka":
        return "results/benchmarking/artic/minion/medaka/{sample}/{sample}.merged.vcf"

    elif pipeline == "artic-nanopolish":
        return (
            "results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.merged.vcf"
        )

    elif pipeline == "ncov2019-artic-nf":
        return "results/benchmarking/ncov2019_artic_nf/illumina/{sample}/ncovIllumina_sequenceAnalysis_callVariants/{sample}.variants.vcf"

    elif pipeline == "ncov2019-artic-nf-medaka":
        return "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{{sample}}_{barcode}.vcf".format(
            barcode=get_barcode(wildcards)
        )

    elif pipeline == "ncov2019-artic-nf-nanopolish":
        return "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{{sample}}-{barcode}/articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish/{{sample}}_{barcode}.merged.vcf".format(
            barcode=get_barcode(wildcards)
        )

    elif pipeline == "nf-core-viralrecon":
        return "results/benchmarking/nf-core-viralrecon/illumina/{sample}/{sample}.vcf"

    elif pipeline == "nf-core-viralrecon-nanopolish":
        return "results/benchmarking/nf-core-viralrecon/nanopore/nanopolish/{sample}/nanopolish/{sample}.merged.vcf"

    elif pipeline == "nf-core-viralrecon-medaka":
        return "results/benchmarking/nf-core-viralrecon/nanopore/medaka/{sample}/medaka/{sample}.merged.vcf"

    elif pipeline == "covpipe":
        return "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/intermediate_data/04_variant_calling/{covpipe_name}/{covpipe_name}.vcf".format(
            covpipe_name=get_covpipe_name_for_sample(wildcards)
        )

    elif pipeline == "havoc":
        return ("results/benchmarking/havoc/{sample}/{sample}.fixed.vcf",)

    elif pipeline == "v-pipe":
        return "results/benchmarking/v-pipe/{sample}/work/samples/{sample}/20200102/variants/SNVs/fixed-vcf/snvs.vcf"

    elif pipeline == "signal":
        return "results/benchmarking/SIGNAL/{sample}/results_dir/{sample}/freebayes/{sample}.variants.norm.vcf"

    elif pipeline == "snakelines":
        return "results/benchmarking/snakelines/{sample}/variant/sars_cov_2-wgs/original/{sample}.vcf"

    elif pipeline == "uncovar":
        return "results/{date}/filtered-calls/ref~main/{{sample}}.subclonal.high+moderate-impact.vcf".format(
            date=get_date_for_sample(wildcards)
        )

    elif pipeline == "sanger":
        return "results/benchmarking/sanger/variant-calls/{sample}.vcf"


def get_sanger_files_for_sample(wildcards):
    return pep.sample_table.loc[wildcards.sample]["sanger"].split(";")
