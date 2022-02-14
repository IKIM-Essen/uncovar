# PIPELINES = {"nanopore": ["uncovar"], "illumina": ["uncovar", "covpipe"]}
PIPELINES = {
    "nanopore": {
        "artic-medaka": {
            "outdir": "results/benchmarking/artic/minion/medaka/{sample}/",
            "vcf": "results/benchmarking/artic/minion/medaka/{sample}/{sample}.merged.vcf",
            "consensus": "results/benchmarking/artic/minion/medaka/{sample}/{sample}.consensus.fasta",
        },
        "artic-nanopolish": {
            "outdir": "results/benchmarking/artic/minion/nanopolish/{sample}/",
            "vcf": "results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.merged.vcf",
            "consensus": "results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.consensus.fasta",
        },
        "ncov2019-artic-nf-medaka": {
            "outdir": "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{{sample}}-{barcode}/",
            "vcf": "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{{sample}}_{barcode}.vcf",
            "consensus": "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{{sample}}-{barcode}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/{{sample}}_{barcode}.consensus.fasta",
        },
        "ncov2019-artic-nf-nanopolish": {
            "outdir": "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{{sample}}-{barcode}/",
            "vcf": "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{{sample}}-{barcode}/articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish/{{sample}}_{barcode}.merged.vcf",
            "consensus": "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{{sample}}-{barcode}/articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish/{{sample}}_{barcode}.consensus.fasta",
        },
        "nf-core-viralrecon-nanopolish": {
            "outdir": "results/benchmarking/nf-core-viralrecon/nanopore/nanopolish/{sample}",
            "vcf": "results/benchmarking/nf-core-viralrecon/nanopore/nanopolish/{sample}/nanopolish/{sample}.merged.vcf",
            "consensus": "results/benchmarking/nf-core-viralrecon/nanopore/nanopolish/{sample}/nanopolish/{sample}.consensus.fasta",
            "pangolin": "results/benchmarking/nf-core-viralrecon/nanopore/nanopolish/{sample}/nanopolish/pangolin/{sample}.pangolin.csv",
        },
        "nf-core-viralrecon-medaka": {
            "outdir": "results/benchmarking/nf-core-viralrecon/nanopore/medaka/{sample}/",
            "vcf": "results/benchmarking/nf-core-viralrecon/nanopore/medaka/{sample}/medaka/{sample}.merged.vcf",
            "consensus": "results/benchmarking/nf-core-viralrecon/nanopore/medaka/{sample}/medaka/{sample}.consensus.fasta",
            "pangolin": "results/benchmarking/nf-core-viralrecon/nanopore/medaka/{sample}/medaka/pangolin/{sample}.pangolin.csv",
        },
        # "poreCov": {
        #     "outdir": "results/benchmarking/poreCov/{sample}/",
        #     "consensus": "results/benchmarking/poreCov/{sample}/2.Genomes/all_consensus_sequences/{sample}.consensus.fasta",
        #     "lineage_call": "results/benchmarking/poreCov/{sample}/3.Lineages_Clades_Mutations/{sample}/lineage_report_{sample}.csv",
        # },
        # "uncovar": {
        #     "outdir": [],
        #     "vcf": "results/{date}/filtered-calls/ref~main/{{sample}}.subclonal.high+moderate-impact.vcf",
        # },
        "sanger": {
            "outdir": [],
            "vcf": "results/benchmarking/sanger/fixed-genotype/{sample}.vcf",
        },
    },
    "illumina": {
        "covpipe": {
            "outdir": "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/",
            "vcf": "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/intermediate_data/04_variant_calling/{covpipe_name}/{covpipe_name}.vcf",
            "consensus": "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/consensuses_masked/{covpipe_name}.masked_consensus.fasta",
            "pangolin": "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/intermediate_data/06_lineages/all_samples.lineage.txt",
        },
        "havoc": {
            "outdir": "results/benchmarking/havoc/{{sample}}/data/{havoc_name}/",
            "vcf": "results/benchmarking/havoc/{sample}/{sample}.fixed.vcf",
            "consensus": "results/benchmarking/havoc/{{sample}}/data/{havoc_name}/{havoc_name}_consensus.fa",
            "pangolin": "results/benchmarking/havoc/{{sample}}/data/{havoc_name}/{havoc_name}_pangolin_lineage.csv",
        },
        "ncov2019-artic-nf": {
            "outdir": "results/benchmarking/ncov2019_artic_nf/illumina/{sample}/",
            "vcf": "results/benchmarking/ncov2019_artic_nf/illumina/{sample}/ncovIllumina_sequenceAnalysis_callVariants/{sample}.variants.vcf",
            "consensus": "results/benchmarking/ncov2019_artic_nf/illumina/{sample}/ncovIllumina_sequenceAnalysis_makeConsensus/{sample}.primertrimmed.consensus.fa",
        },
        "nf-core-viralrecon": {
            "outdir": "results/benchmarking/nf-core-viralrecon/illumina/{sample}",
            "vcf": "results/benchmarking/nf-core-viralrecon/illumina/{sample}.vcf",
            "consensus": "results/benchmarking/nf-core-viralrecon/illumina/{sample}/variants/bcftools/consensus/{sample}.consensus.fa",
            "pangolin": "results/benchmarking/nf-core-viralrecon/illumina/{sample}/variants/bcftools/pangolin/{sample}.pangolin.csv",
            "de_novo_assembly": "results/benchmarking/nf-core-viralrecon/illumina/{sample}/assembly/spades/rnaviral/{sample}.contigs.fa",
        },
        "signal": {
            "outdir": "results/benchmarking/SIGNAL/{sample}/results_dir",
            "vcf": "results/benchmarking/SIGNAL/{sample}/results_dir/{sample}/freebayes/{sample}.variants.norm.vcf",
            "consensus": "results/benchmarking/SIGNAL/{sample}/results_dir/all_freebayes_genomes.fa",
            "pangolin": "results/benchmarking/SIGNAL/{sample}/results_dir/freebayes_lineage_assignments.tsv",
        },
        "snakelines": {
            "outdir": "results/benchmarking/snakelines/{sample}/",
            "vcf": "results/benchmarking/snakelines/{sample}/variant/sars_cov_2-wgs/original/{sample}.vcf",
            "consensus": "results/benchmarking/snakelines/{sample}/report/public/01-example/{sample}/consensus-sars_cov_2-wgs.fa",
            "pangolin": "results/benchmarking/snakelines/{sample}/report/public/01-example/{sample}/lineage_report-sars_cov_2-wgs.csv",
        },
        # "uncovar": {
        #     "outdir": [],
        #     "vcf": "results/{date}/filtered-calls/ref~main/{{sample}}.subclonal.high+moderate-impact.vcf",
        # },
        "v-pipe": {
            "outdir": "results/benchmarking/v-pipe/{sample}/work",
            "vcf": "results/benchmarking/v-pipe/{sample}/work/samples/{sample}/20200102/variants/SNVs/snvs.vcf",
            "consensus": "results/benchmarking/v-pipe/{sample}/work/samples/{sample}/20200102/references/ref_majority.fasta",
        },
        "sanger": {
            "outdir": [],
            "vcf": "results/benchmarking/sanger/fixed-genotype/{sample}.vcf",
        },
    },
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


def get_output_from_pipline(key):
    def inner(wildcards):
        try:
            path = PIPELINES["illumina"][wildcards.workflow][key]
        except KeyError:
            path = PIPELINES["nanopore"][wildcards.workflow][key]

        if "{barcode}" in path:
            path = path.format(barcode=get_barcode(wildcards))

        if "{covpipe_name}" in path:
            path = path.format(covpipe_name=get_covpipe_name_for_sample(wildcards))

        if "{date}" in path:
            path = path.format(date=get_date_for_sample(wildcards))
        if "{havoc_name}" in path:
            path = path.format(havoc_name=wildcards.sample.split("_")[0])
        return path

    return inner


def get_sanger_files_for_sample(wildcards):
    return pep.sample_table.loc[wildcards.sample]["sanger"].split(";")


def get_benchmark_paths_by_tech(path, tech, samples):
    return expand(
        path,
        workflow=PIPELINES[tech].keys(),
        sample=samples,
    )


def get_benchmark_path(path):
    def inner(wildcards):
        return get_benchmark_paths_by_tech(
            path, "nanopore", get_nanopore_samples(wildcards)
        ) + get_benchmark_paths_by_tech(
            path, "illumina", get_illumina_samples(wildcards)
        )

    return inner


def get_benchmark_platforms(wildcards):
    return ["nanopore"] * len(
        get_benchmark_paths_by_tech(
            "{workflow},{sample}", "nanopore", get_nanopore_samples(wildcards)
        )
    ) + ["illumina"] * len(
        get_benchmark_paths_by_tech(
            "{workflow},{sample}", "illumina", get_illumina_samples(wildcards)
        )
    )
