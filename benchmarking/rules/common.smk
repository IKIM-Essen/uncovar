PIPELINES = {
    "nanopore": {
        "artic-medaka": {
            "outdir": "results/benchmarking/artic/minion/medaka/{sample}/",
            "vcf": "results/benchmarking/artic/minion/medaka/{sample}/{sample}.merged.vcf",
            "consensus": "results/benchmarking/artic/minion/medaka/{sample}/{sample}.consensus.fasta",
            "time": "benchmarks/artic_medaka/{sample}.all.benchmark.txt",
        },
        "artic-nanopolish": {
            "outdir": "results/benchmarking/artic/minion/nanopolish/{sample}/",
            "vcf": "results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.merged.vcf",
            "consensus": "results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.consensus.fasta",
            "time": "benchmarks/artic_nanopolish/{sample}.all.benchmark.txt",
        },
        # "ncov2019-artic-nf-medaka": { # Must specify --medaka-model if using the --medaka workflow. But this is not possible as ncov2019 does not provide an interface for tht.
        #     "outdir": "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{{sample}}-{barcode}/",
        #     "vcf": "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{{sample}}_{barcode}.vcf",
        #     "consensus": "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{{sample}}-{barcode}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/{{sample}}_{barcode}.consensus.fasta",
        #     "time": "benchmarks/ncov2019_artic_nf_medaka/{{sample}}~{barcode}.benchmark.txt"
        # },
        "ncov2019-artic-nf-nanopolish": {
            "outdir": "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{{sample}}-{barcode}/",
            "vcf": "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{{sample}}-{barcode}/articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish/{{sample}}_{barcode}.merged.vcf",
            "consensus": "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{{sample}}-{barcode}/articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish/{{sample}}_{barcode}.consensus.fasta",
            "time": "benchmarks/ncov2019_artic_nf_nanopolish/{{sample}}~{barcode}.benchmark.txt",
        },
        "nf-core-viralrecon-nanopolish": {
            "outdir": "results/benchmarking/nf-core-viralrecon/nanopore/nanopolish/{sample}",
            "vcf": "results/benchmarking/nf-core-viralrecon/nanopore/nanopolish/{sample}/nanopolish/{sample}.merged.vcf",
            "consensus": "results/benchmarking/nf-core-viralrecon/nanopore/nanopolish/{sample}/nanopolish/{sample}.consensus.fasta",
            "pangolin": "results/benchmarking/nf-core-viralrecon/nanopore/nanopolish/{sample}/nanopolish/pangolin/{sample}.pangolin.csv",
            "time": "benchmarks/nf_core_viralrecon_nanopolish/{sample}.benchmark.txt",
        },
        "nf-core-viralrecon-medaka": {
            "outdir": "results/benchmarking/nf-core-viralrecon/nanopore/medaka/{sample}/",
            "vcf": "results/benchmarking/nf-core-viralrecon/nanopore/medaka/{sample}/medaka/{sample}.merged.vcf",
            "consensus": "results/benchmarking/nf-core-viralrecon/nanopore/medaka/{sample}/medaka/{sample}.consensus.fasta",
            "pangolin": "results/benchmarking/nf-core-viralrecon/nanopore/medaka/{sample}/medaka/pangolin/{sample}.pangolin.csv",
            "time": "benchmarks/nf_core_viralrecon_medaka/{sample}.benchmark.txt",
        },
        "poreCov": {  # no vcf output
            "outdir": "results/benchmarking/poreCov/{sample}/",
            "consensus": "results/benchmarking/poreCov/{sample}/2.Genomes/all_consensus_sequences/{sample}.consensus.fasta",
            "lineage_call": "results/benchmarking/poreCov/{sample}/3.Lineages_Clades_Mutations/{sample}/lineage_report_{sample}.csv",
            "time": "benchmarks/porecov/{sample}.benchmark.txt",
        },
        "uncovar": {
            "outdir": [],
            "vcf": "results/benchmarking/UnCoVar/benchmark-result/{sample}.vcf",
            "time": "benchmarks/uncovar/{{sample}}~{date}~consensus.benchmark.txt",
            "pangolin": "results/benchmarking/UnCoVar/results/{date}/tables/strain-calls/{{sample}}.consensus.strains.pangolin.csv",
        },
        "sanger": {
            "outdir": [],
            "vcf": "results/benchmarking/sanger/fixed-genotype/{sample}.vcf",
        },
    },
    "illumina": {
        "covpipe": {
            "outdir": "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/",
            "vcf": "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/intermediate_data/04_variant_calling/{{sample}}_/{{sample}}_.vcf",
            "consensus": "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/consensuses_masked/{{sample}}_.masked_consensus.fasta",
            "pangolin": "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/intermediate_data/06_lineages/{{sample}}_/{{sample}}_.lineage.txt",
            "time": "benchmarks/covpipe/{{sample}}~{covpipe_name}.benchmark.txt",
        },
        "havoc": {
            "outdir": "results/benchmarking/havoc/{{sample}}/data/{havoc_name}/",
            "vcf": "results/benchmarking/havoc/{sample}/{sample}.fixed.vcf",
            "consensus": "results/benchmarking/havoc/{{sample}}/data/{havoc_name}/{havoc_name}_consensus.fa",
            "pangolin": "results/benchmarking/havoc/{{sample}}/data/{havoc_name}/{havoc_name}_pangolin_lineage.csv",
            "time": "benchmarks/havoc/{{sample}}~{havoc_name}.benchmark.txt",
        },
        "ncov2019-artic-nf": {
            "outdir": "results/benchmarking/ncov2019_artic_nf/illumina/{sample}/",
            "vcf": "results/benchmarking/ncov2019_artic_nf/illumina/{sample}/ncovIllumina_sequenceAnalysis_callVariants/{sample}.variants.vcf",
            "consensus": "results/benchmarking/ncov2019_artic_nf/illumina/{sample}/ncovIllumina_sequenceAnalysis_makeConsensus/{sample}.primertrimmed.consensus.fa",
            "time": "benchmarks/ncov2019_artic_nf_illumina/{sample}.benchmark.txt",
        },
        "nf-core-viralrecon": {
            "outdir": "results/benchmarking/nf-core-viralrecon/illumina/{sample}",
            "vcf": "results/benchmarking/nf-core-viralrecon/illumina/{sample}.vcf",
            "consensus": "results/benchmarking/nf-core-viralrecon/illumina/{sample}/variants/bcftools/consensus/{sample}.consensus.fa",
            "pangolin": "results/benchmarking/nf-core-viralrecon/illumina/{sample}/variants/bcftools/pangolin/{sample}.pangolin.csv",
            "de_novo_assembly": "results/benchmarking/nf-core-viralrecon/illumina/{sample}/assembly/spades/rnaviral/{sample}.contigs.fa",
            "time": "benchmarks/nf_core_viralrecon_illumina/{sample}.benchmark.txt",
        },
        # "signal": { # failes bc of pangolin update, which uses GiTHub API, which is rate limitet to 60 requests
        #     "outdir": "results/benchmarking/SIGNAL/{sample}/results_dir",
        #     "vcf": "results/benchmarking/SIGNAL/{sample}/results_dir/{sample}/freebayes/{sample}.variants.norm.vcf",
        #     "consensus": "results/benchmarking/SIGNAL/{sample}/results_dir/all_freebayes_genomes.fa",
        #     "pangolin": "results/benchmarking/SIGNAL/{sample}/results_dir/freebayes_lineage_assignments.tsv",
        #     "time":"benchmarks/signal/{sample}.benchmark.txt"
        # },
        "snakelines": {
            "outdir": "results/benchmarking/snakelines/{sample}/",
            "vcf": "results/benchmarking/snakelines/{sample}/variant/sars_cov_2-wgs/original/{sample}.vcf",
            "consensus": "results/benchmarking/snakelines/{sample}/report/public/01-example/{sample}/consensus-sars_cov_2-wgs.fa",
            "pangolin": "results/benchmarking/snakelines/{sample}/report/public/01-example/{sample}/lineage_report-sars_cov_2-wgs.csv",
            "time": "benchmarks/snakeline/{sample}.benchmark.txt",
        },
        "uncovar": {
            "outdir": [],
            "vcf": "results/benchmarking/UnCoVar/benchmark-result/{sample}.vcf",
            "time": "benchmarks/uncovar/{{sample}}~{date}~polished.benchmark.txt",
            "pangolin": "results/benchmarking/UnCoVar/results/{date}/tables/strain-calls/{{sample}}.polished.strains.pangolin.csv",
        },
        "v-pipe": {
            "outdir": [
                "results/benchmarking/v-pipe/{sample}/samples/{sample}/20200102/variants",
                "results/benchmarking/v-pipe/{sample}/samples/{sample}/20200102/references",
            ],
            "vcf": "results/benchmarking/v-pipe/fixed-vcf/{sample}.vcf",
            "consensus": "results/benchmarking/v-pipe/{sample}/samples/{sample}/20200102/references/ref_majority.fasta",
            "time": "benchmarks/v_pipe/{sample}.benchmark.txt",
        },
        "sanger": {
            "outdir": [],
            "vcf": "results/benchmarking/sanger/fixed-genotype/{sample}.vcf",
        },
    },
}

ILLUMINA = "illumina"
ONT = "ont"

COVERAGES = {
    "none": 0,
    "low": 1,
    "high": 20,
}


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_output_dir(wildcards, output):
    return os.path.dirname(output[0])


def get_technology(wildcards, sample=None):
    if sample is None:
        sample = wildcards.sample

    return pep.sample_table.loc[sample]["technology"]


def is_illumina(wildcards, sample=None):
    if sample is None:
        return get_technology(wildcards) == ILLUMINA
    return get_technology(None, sample) == ILLUMINA


def is_ont(wildcards, sample=None):
    if sample is None:
        return get_technology(wildcards) == ONT
    return get_technology(None, sample) == ONT


def get_fastqs(wildcards):
    if is_illumina(wildcards):
        return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]
    elif is_ont(wildcards):
        return pep.sample_table.loc[wildcards.sample][["fq1"]]


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


def check_path_for_workflow_wildcard(path, wildcards):
    if "{barcode}" in path:
        path = path.format(barcode=get_barcode(wildcards))

    if "{covpipe_name}" in path:
        path = path.format(covpipe_name=get_covpipe_name_for_sample(wildcards))

    if "{date}" in path:
        path = path.format(date=get_date_for_sample(wildcards))

    if "{havoc_name}" in path:
        path = path.format(havoc_name=wildcards.sample.split("_")[0])

    return path


def get_output_from_pipline(key):
    def inner(wildcards):
        try:
            path = PIPELINES["illumina"][wildcards.workflow][key]
        except KeyError:
            path = PIPELINES["nanopore"][wildcards.workflow][key]

        return check_path_for_workflow_wildcard(path, wildcards)

    return inner


def get_sanger_files_for_sample(wildcards):
    return pep.sample_table.loc[wildcards.sample]["sanger"].split(";")


def get_benchmark_paths_by_tech(
    path, tech, samples, remove=None, pop_samples=None, key=None
):

    workflows = list(PIPELINES[tech].keys())

    if key is not None:
        for workflow in workflows:
            if key not in PIPELINES[tech][workflow].keys():
                workflows.remove(workflow)

    if pop_samples is not None:
        samples = samples.tolist()
        [samples.remove(sample) for sample in pop_samples if sample in samples]

    if remove in workflows:
        workflows.remove(remove)

    return expand(
        path,
        workflow=workflows,
        sample=samples,
    )


def get_benchmark_path(path, remove=None, pop_sample=None):
    def inner_get_benchmark_path(wildcards):
        return get_benchmark_paths_by_tech(
            path, "nanopore", get_nanopore_samples(wildcards), remove, pop_sample
        ) + get_benchmark_paths_by_tech(
            path, "illumina", get_illumina_samples(wildcards), remove, pop_sample
        )

    return inner_get_benchmark_path


def get_benchmark_platforms(remove=None, pop_sample=None):
    def inner(wildcards):
        return ["nanopore"] * len(
            get_benchmark_paths_by_tech(
                "{workflow},{sample}",
                "nanopore",
                get_nanopore_samples(wildcards),
                remove,
                pop_sample,
            )
        ) + ["illumina"] * len(
            get_benchmark_paths_by_tech(
                "{workflow},{sample}",
                "illumina",
                get_illumina_samples(wildcards),
                remove,
                pop_sample,
            )
        )

    return inner


def get_workflow_output(wildcards):
    path = PIPELINES[wildcards.tech][wildcards.workflow][wildcards.key]
    return check_path_for_workflow_wildcard(path, wildcards)


def get_all_outputs(wildcards):
    path = (
        "results/benchmarking/backups/{key}/{tech}/{workflow}/{sample}.some.extension"
    )

    paths = []
    for tech, tech_dict in PIPELINES.items():
        for workflow, workflow_dict in tech_dict.items():
            for key, _dict_path in workflow_dict.items():
                if key != "outdir":
                    if tech == "illumina":
                        samples = get_illumina_samples(wildcards)
                    if tech == "nanopore":
                        samples = get_nanopore_samples(wildcards)

                    for sample in samples:
                        paths.append(
                            path.format(
                                key=key, tech=tech, workflow=workflow, sample=sample
                            )
                        )
    return paths


def get_happy_output(path):
    def inner(wildcards):
        with checkpoints.get_samples_with_multiallelic_calls.get(
            source=wildcards.source, cov=wildcards.cov
        ).output[0].open() as f:
            multiallelic_calls = [line.split("\t")[0] for line in f.readlines()]

        with checkpoints.extract_truth_without_calls.get(
            source=wildcards.source, cov=wildcards.cov
        ).output[0].open() as f:
            truth_without_calls = [line.strip() for line in f.readlines()]

        samples_to_exclude = list(set(multiallelic_calls + truth_without_calls))

        return get_benchmark_paths_by_tech(
            path=path,
            tech="nanopore",
            samples=get_nanopore_samples(wildcards),
            remove="sanger",
            pop_samples=samples_to_exclude,
            key="vcf",
        ) + get_benchmark_paths_by_tech(
            path=path,
            tech="illumina",
            samples=get_illumina_samples(wildcards),
            remove="sanger",
            pop_samples=samples_to_exclude,
            key="vcf",
        )

    return inner


def get_happy_platform_data(remove):
    def inner(wildcards):
        with checkpoints.get_samples_with_multiallelic_calls.get(
            source=wildcards.source, cov=wildcards.cov
        ).output[0].open() as f:
            multiallelic_calls = [line.split("\t")[0] for line in f.readlines()]

        with checkpoints.extract_truth_without_calls.get(
            source=wildcards.source, cov=wildcards.cov
        ).output[0].open() as f:
            truth_without_calls = [line.strip() for line in f.readlines()]

        samples_to_exclude = list(set(multiallelic_calls + truth_without_calls))

        return ["nanopore"] * len(
            get_benchmark_paths_by_tech(
                "{workflow},{sample}",
                "nanopore",
                get_nanopore_samples(wildcards),
                remove,
                samples_to_exclude,
            )
        ) + ["illumina"] * len(
            get_benchmark_paths_by_tech(
                "{workflow},{sample}",
                "illumina",
                get_illumina_samples(wildcards),
                remove,
                samples_to_exclude,
            )
        )

    return inner


def get_output_of_pipelines(path, output):
    def inner(wildcards):
        expanded_paths = []
        for tech, workflow_dict in PIPELINES.items():
            for workflow, output_dict in workflow_dict.items():
                if output in output_dict.keys():
                    if tech == "illumina":
                        samples = get_illumina_samples(wildcards)
                    elif tech == "nanopore":
                        samples = get_nanopore_samples(wildcards)

                    paths = expand(
                        path,
                        output_type=output,
                        tech=tech,
                        workflow=workflow,
                        sample=samples,
                    )
                    expanded_paths.append(paths)

        # flatten list of lists
        expanded_paths = [item for sublist in expanded_paths for item in sublist]

        assert len(expanded_paths) > 0, f"No files found for {output}"

        return expanded_paths

    return inner


def get_paths_of_input_files(wildcards):
    if is_illumina(wildcards):
        return get_fastqs(wildcards)[0]
    elif is_ont(wildcards):
        fastq_path = get_fastq_pass_path_barcode(wildcards)
        return [
            os.path.join(fastq_path, f)
            for f in os.listdir(fastq_path)
            if os.path.isfile(os.path.join(fastq_path, f))
        ]

    raise NotImplementedError(f"No file found for {wildcards}")


def get_mosdepth_quantize():
    return ":".join(map(str, sorted(COVERAGES.values()))) + ":"


def get_io_prefix(getter):
    def inner(wildcards, input, output):

        return getter(input, output).split(".")[0]

    return inner


def get_cov_label(wildcards):
    lower, upper = get_cov_interval(wildcards.cov)
    if upper:
        return f"{lower}:{upper}"
    return f"{lower}:inf"


def get_cov_interval(name):
    threshold = COVERAGES[name]
    upper_bound = None

    greater = [cov for cov in COVERAGES.values() if cov > threshold]
    if greater:
        upper_bound = min(greater)

    return threshold, upper_bound


def get_adapters(wildcards):
    try:
        samples = pep.sample_table
        samples.dropna(subset=["adapters"], inplace=True)
        adapters = samples.loc[wildcards.sample]["adapters"]

        # predefined adapters
        # https://lifesciences.tecan.com/revelo-rna-seq-library-prep-kit
        if adapters == "revelo-rna-seq":
            return "--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
        # https://www.nimagen.com/shop/products/rc-cov096/easyseq-sars-cov-2-novel-coronavirus-whole-genome-sequencing-kit
        elif adapters == "nimagen-easy-seq":
            return "--adapter_sequence GCGAATTTCGACGATCGTTGCATTAACTCGCGAA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

    # If there is no column namend "adapters" in the sample sheet
    # return the adapters defined in the config file
    except KeyError:
        return config["kit-adapters"]


wildcard_constraints:
    sample="[^/.]+",
