# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

from pathlib import Path
import pandas as pd
import re
import random
from snakemake.utils import validate


VARTYPES = ["SNV", "MNV", "INS", "DEL", "REP", "INV", "DUP"]
ILLUMINA_VARRANGE = ["small", "structural"]
ONT_VARRANGE = ["homopolymer-medaka", "homopolymer-longshot"]
ION_VARRANGE = ["small", "structural"]

# clear text / content of flag "technology" in sample sheet
ILLUMINA = "illumina"
ONT = "ont"
ION_TORRENT = "ion"

# for benchmarking rules
BENCHMARK_PREFIX = "benchmark-sample-"
NON_COV2_TEST_PREFIX = "non-cov2-"
MIXTURE_PREFIX = "mixture-sample-"
MIXTURE_PART_INDICATOR = "_MIX_"
MIXTURE_PERCENTAGE_INDICATOR = "_PERC_"
BENCHMARK_DATE_WILDCARD = "benchmarking"
READ_TEST_PREFIX = "read-sample-"
READ_NUMBER_INDICATOR = "_READ_NUMBER_"
READ_LENGTH_INDICATOR = "_READ_LENGTH_"
READ_STATE_INDICATOR = "_STATE_"


configfile: "config/config.yaml"


validate(config, "../schemas/config.schema.yaml")

validate(pep.sample_table, "../schemas/samples.schema.yaml")


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_dates():
    return list(pep.sample_table["date"].values)


def get_samples_for_date(date, filtered=False):
    # select samples with given date
    df = pep.sample_table
    df = df[df["date"] == date]

    samples_of_run = list(df["sample_name"].values)

    # filter
    if filtered:
        with checkpoints.quality_filter.get(date=date).output.passed_filter.open() as f:

            passend_samples = []
            for line in f:
                passend_samples.append(line.strip())

        filtered_samples = [
            sample for sample in samples_of_run if sample in passend_samples
        ]

        if not filtered_samples:
            raise ValueError(
                "List of filtered samples is empty. Perhaps no samples of run {} passed the quality criteria.".format(
                    date
                )
            )

        return filtered_samples

    # unfiltered
    else:
        return samples_of_run


def get_all_run_dates():
    sorted_list = list(pep.sample_table["date"].unique())
    sorted_list.sort()
    return sorted_list


def get_latest_run_date():
    return pep.sample_table["date"].max()


def get_samples_before_date(wildcards):
    return list(
        pep.sample_table[pep.sample_table["date"] <= wildcards.date][
            "sample_name"
        ].values
    )


def get_dates_before_date(wildcards):
    return list(
        pep.sample_table[pep.sample_table["date"] <= wildcards.date]["date"].values
    )


def get_technology(wildcards, sample=None):
    benchmark_technology = ILLUMINA

    if sample is None:
        sample = wildcards.sample

    if is_benchmark_data(sample):
        return benchmark_technology
    return pep.sample_table.loc[sample]["technology"]


def is_ont(wildcards, sample=None):
    if sample is None:
        return get_technology(wildcards) == ONT
    return get_technology(None, sample) == ONT


def is_illumina(wildcards, sample=None):
    if sample is None:
        return get_technology(wildcards) == ILLUMINA
    return get_technology(None, sample) == ILLUMINA


def is_ion_torrent(wildcards, sample=None):
    """Returns if the sample was sequenced with the Ion Torrent platform."""
    if sample is None:
        return get_technology(wildcards) == ION_TORRENT
    return get_technology(None, sample) == ION_TORRENT


def has_pseudo_assembly(wildcards, sample=None):
    """Returns if a pseudo-assembly should be created for the sample."""
    if sample is None:
        return is_illumina(wildcards) or is_ion_torrent(wildcards)
    return is_illumina(None, sample) or is_ion_torrent(None, sample)


def has_consensus_assembly(wildcards, sample=None):
    """Returns if a consensus-assembly should be created for the sample."""
    if sample is None:
        return is_ont(wildcards)
    return is_ont(None, sample)


def is_single_end(wildcards, sample=None):
    """Returns if the sample was sequenced with single end technology."""
    if sample is None:
        return is_ont(wildcards) or is_ion_torrent(wildcards)
    return is_ont(None, sample) or is_ion_torrent(None, sample)


def get_fastqs(wildcards):
    if wildcards.sample.startswith(BENCHMARK_PREFIX):
        # this is a simulated benchmark sample, do not look up FASTQs in the sample sheet
        accession = wildcards.sample[len(BENCHMARK_PREFIX) :]
        return expand(
            "resources/benchmarking/{accession}/reads.{read}.fastq.gz",
            accession=accession,
            read=[1, 2],
        )
    # non-sars-cov2-genome test
    if wildcards.sample.startswith(NON_COV2_TEST_PREFIX):
        # this is for testing non-sars-cov2-genomes
        accession = wildcards.sample[len(NON_COV2_TEST_PREFIX) :]
        return expand(
            "resources/benchmarking/{accession}/reads.{read}.fastq.gz",
            accession=accession,
            read=[1, 2],
        )
    # mixture
    if wildcards.sample.startswith(MIXTURE_PREFIX):
        mixture = wildcards.sample[len(MIXTURE_PREFIX) :]
        return expand(
            "resources/mixtures/{mixtures}/reads.{read}.fastq.gz",
            mixtures=mixture,
            read=[1, 2],
        )
    # read benchmark
    if wildcards.sample.startswith(READ_TEST_PREFIX):
        return expand(
            "resources/benchmarking/{accession}/reads.{read}.fastq.gz",
            accession=wildcards.sample,
            read=[1, 2],
        )

    # default case, look up FASTQs in the sample sheet
    if is_illumina(wildcards):
        return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]
    elif is_ont(wildcards):
        return pep.sample_table.loc[wildcards.sample][["fq1"]]
    elif is_ion_torrent(wildcards):
        return pep.sample_table.loc[wildcards.sample][["fq1"]]


def get_resource(name):
    return str((Path(workflow.snakefile).parent.parent.parent / "resources") / name)


def get_report_input(pattern):
    def inner(wildcards):
        return expand(pattern, sample=get_report_samples(wildcards), **wildcards)

    return inner


def get_report_bcfs(wildcards, input):
    """Return paths to BCF files for reporting."""
    return expand(
        "{sample}={bcf}", zip, sample=get_report_samples(wildcards), bcf=input.bcfs
    )


def get_report_bams(wildcards, input):
    return expand(
        "{sample}:{sample}={bam}",
        zip,
        sample=get_report_samples(wildcards),
        bam=input.bams,
    )


def get_report_samples(wildcards):
    return (
        get_samples_for_date(wildcards.date)
        if wildcards.target == "all"
        else [wildcards.target]
    )


def get_merge_calls_input(suffix):
    def inner(wildcards):
        return expand(
            "results/{{date}}/filtered-calls/ref~{{reference}}/{{sample}}.{{clonality}}.{{filter}}.{vartype}.fdr-controlled{suffix}",
            suffix=suffix,
            vartype=VARTYPES,
        )

    return inner


def get_strain_accessions(wildcards):
    with checkpoints.get_strain_accessions.get().output[0].open() as f:
        # Get genomes for benchmarking from config
        accessions = config.get("testing", {}).get("benchmark-genomes", [])
        if not accessions:
            accessions = pd.read_csv(f, squeeze=True)
        try:
            accessions = accessions[: config["testing"]["limit-strain-genomes"]]
        except KeyError:
            # take all strain genomes
            pass
        return accessions


def get_non_cov2_accessions():
    accessions = config.get("non_cov2_genomes", [])
    return accessions


def load_strain_genomes(f):
    strain_genomes = pd.read_csv(f, squeeze=True).to_list()
    strain_genomes.append("resources/genomes/main.fasta")
    return expand("{strains}", strains=strain_genomes)


def get_strain_genomes(wildcards):
    # Case 1: take custom genomes from gisaid
    if not config.get("testing", {}).get("use-genbank", False):
        if config["strain-calling"]["use-gisaid"]:
            # use genomes extracted from gisaid provision
            with checkpoints.extract_strain_genomes_from_gisaid.get(
                date=wildcards.date
            ).output[0].open() as f:
                return load_strain_genomes(f)
        # use genomes from genbank
        with checkpoints.get_lineages_for_non_gisaid_based_calling.get(
            date=wildcards.date
        ).output[0].open() as f:
            return load_strain_genomes(f)

    # Case 2: for benchmarking (no strain-calling/genomes in config file)
    # take genomes from genbank
    accessions = get_strain_accessions(wildcards)
    return expand("resources/genomes/{accession}.fasta", accession=accessions)


def get_strain_signatures(wildcards):
    return expand(
        "resources/genomes/{accession}.sig", accession=get_strain_accessions(wildcards)
    )


def get_benchmark_results(wildcards):
    accessions = get_strain_accessions(wildcards)
    return expand(
        "results/benchmarking/tables/strain-calls/benchmark-sample-{accession}.strains.kallisto.tsv",
        accession=accessions,
    )


def get_assembly_comparisons(bams=True):
    def inner(wildcards):
        accessions = get_strain_accessions(wildcards)
        pattern = (
            "results/benchmarking/assembly/{{assembly_type}}/{accession}.bam"
            if bams
            else "resources/genomes/{accession}.fasta"
        )
        return expand(
            pattern,
            accession=accessions,
        )

    return inner


def get_assembly_result(wildcards):
    if wildcards.assembly_type == "assembly":
        return (
            "results/benchmarking/contigs/polished/benchmark-sample-{accession}.fasta"
        )
    elif wildcards.assembly_type == "pseudoassembly":
        return "results/benchmarking/contigs/pseudoassembled/benchmark-sample-{accession}.fasta"
    else:
        raise ValueError(
            f"unexpected value for wildcard assembly_type: {wildcards.assembly_type}"
        )


def get_non_cov2_calls(from_caller="pangolin"):
    accessions = get_non_cov2_accessions()
    pattern = (
        "results/benchmarking/tables/strain-calls/non-cov2-{accession}.polished.strains.pangolin.csv"
        if from_caller == "pangolin"
        else "results/benchmarking/tables/strain-calls/non-cov2-{accession}.strains.kallisto.tsv"
        if from_caller == "kallisto"
        else []
    )

    if not pattern:
        raise NameError(f"Caller {from_caller} not recognized")

    return expand(pattern, accession=accessions)


def get_reference(suffix=""):
    def inner(wildcards):
        if wildcards.reference == "main":
            # return covid reference genome
            return "resources/genomes/main.fasta{suffix}".format(suffix=suffix)
        elif wildcards.reference == "human":
            # return human reference genome
            return "resources/genomes/human-genome.fna.gz"
        elif wildcards.reference == "main+human":
            return "resources/genomes/main-and-human-genome.fna.gz"
        elif wildcards.reference.startswith("polished-"):
            # return polished contigs
            return "results/{date}/contigs/polished/{sample}.fasta".format(
                sample=wildcards.reference.replace("polished-", ""), **wildcards
            )
        elif wildcards.reference.startswith("consensus-"):
            # return consensus contigs
            return "results/{date}/contigs/consensus/{sample}.fasta".format(
                sample=wildcards.reference.replace("consensus-", ""), **wildcards
            )
        elif wildcards.reference == config["preprocessing"]["amplicon-reference"]:
            # return reference genome of amplicon primers
            return "resources/genomes/{reference}.fasta{suffix}".format(
                reference=config["preprocessing"]["amplicon-reference"], suffix=suffix
            )
        else:
            # return assembly result
            return "results/{date}/contigs/ordered/{reference}.fasta{suffix}".format(
                suffix=suffix, **wildcards
            )

    return inner


def get_bwa_index_prefix(index_paths):
    return os.path.splitext(index_paths[0])[0]


def get_reads(wildcards):
    # alignment against the human reference genome is done with trimmed reads,
    # since this alignment is used to generate the ordered, non human reads
    if (
        wildcards.reference == "human"
        or wildcards.reference == "main+human"
        or wildcards.reference.startswith("polished-")
        or wildcards.reference.startswith("consensus-")
    ):

        illumina_pattern = expand(
            "results/{date}/trimmed/fastp-pe/{sample}.{read}.fastq.gz",
            read=[1, 2],
            **wildcards,
        )

        ont_pattern = expand(
            "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
            **wildcards,
        )

        ion_torrent_pattern = expand(
            "results/{date}/trimmed/fastp-se/{sample}.fastq.gz",
            **wildcards,
        )

        return get_pattern_by_technology(
            wildcards,
            illumina_pattern=illumina_pattern,
            ont_pattern=ont_pattern,
            ion_torrent_pattern=ion_torrent_pattern,
        )

    # theses reads are used to generate the bam file for the BAMclipper and the coverage plot of the main reference
    elif wildcards.reference == config["preprocessing"]["amplicon-reference"]:
        return get_non_human_reads(wildcards)

    # aligments to the references main reference genome,
    # are done with reads, which have undergone the quality control process
    else:
        return get_reads_after_qc(wildcards)


def get_non_human_reads(wildcards):
    return get_pattern_by_technology(
        wildcards,
        illumina_pattern=expand(
            "results/{{date}}/nonhuman-reads/pe/{{sample}}.{read}.fastq.gz",
            read=[1, 2],
        ),
        ont_pattern="results/{date}/nonhuman-reads/se/{sample}.fastq.gz",
        ion_torrent_pattern="results/{date}/nonhuman-reads/se/{sample}.fastq.gz",
    )


def get_reads_after_qc(wildcards, read="both"):
    # in generall: trimmed -> non-human -> clipped reads
    # for shotgun data use non-human reads

    pattern = []

    if is_amplicon_data(wildcards.sample):
        illumina_pattern = expand(
            "results/{date}/read-clipping/fastq/pe/{sample}.{read}.fastq.gz",
            read=[1, 2],
            **wildcards,
        )
        ont_pattern = expand(
            "results/{date}/nonhuman-reads/se/{sample}.fastq.gz", **wildcards
        )
        ion_torrent_pattern = expand(
            "results/{date}/read-clipping/fastq/se/{sample}.fastq", **wildcards
        )

        pattern = get_pattern_by_technology(
            wildcards,
            illumina_pattern=illumina_pattern,
            ont_pattern=ont_pattern,
            ion_torrent_pattern=ion_torrent_pattern,
        )
    # shotgun reads
    else:
        illumina_pattern = expand(
            "results/{date}/nonhuman-reads/pe/{sample}.{read}.fastq.gz",
            read=[1, 2],
            **wildcards,
        )

        ont_pattern = expand(
            "results/{date}/nonhuman-reads/se/{sample}.fastq.gz", **wildcards
        )

        ion_torrent_pattern = expand(
            "results/{date}/nonhuman-reads/se/{sample}.fastq.gz", **wildcards
        )

        pattern = get_pattern_by_technology(
            wildcards,
            illumina_pattern=illumina_pattern,
            ont_pattern=ont_pattern,
            ion_torrent_pattern=ion_torrent_pattern,
        )

    if not pattern:
        raise NotImplementedError(
            f"UnCoVer currently does not support non-amplicon processing for sample {wildcards.sample}"
        )

    if read == "1":
        return pattern[0]
    if read == "2":
        return pattern[1]

    return pattern


def get_min_coverage(wildcards):
    conf = config["quality-criteria"]
    if is_amplicon_data(wildcards.sample):
        return conf["min-depth-with-PCR-duplicates"]
    else:
        return conf["min-depth-without-PCR-duplicates"]


def return_assembler(sample):
    pattern = []
    if is_amplicon_data(sample):
        pattern = get_pattern_by_technology(
            None,
            sample=sample,
            illumina_pattern="{assembler}-pe".format(
                assembler=config["assembly"]["illumina"]["amplicon"]
            ),
            ont_pattern="{assembler}-se".format(
                assembler=config["assembly"]["oxford nanopore"]["amplicon"]
            ),
            ion_torrent_pattern="{assembler}-se".format(
                assembler=config["assembly"]["ion torrent"]["amplicon"]
            ),
        )
    else:
        pattern = get_pattern_by_technology(
            None,
            sample=sample,
            illumina_pattern="{assembler}-pe".format(
                assembler=config["assembly"]["illumina"]["shotgun"]
            ),
            ont_pattern="{assembler}-se".format(
                assembler=config["assembly"]["oxford nanopore"]["shotgun"]
            ),
            ion_torrent_pattern="{assembler}-se".format(
                assembler=config["assembly"]["ion torrent"]["shotgun"]
            ),
        )

    if pattern:
        return pattern

    raise NotImplementedError(
        'No assembler found for technology "{technology}" (sample {sample}).'.format(
            technology=get_technology(wildcards), sample=wildcards.sample
        )
    )


def get_contigs(wildcards, opt_sample=None):
    if "sample" in wildcards.keys():
        return "results/{{date}}/assembly/{{sample}}/{assembler}/{{sample}}.contigs.fasta".format(
            assembler=return_assembler(wildcards.sample)
        )

    # wildcards is only sample name
    return (
        "results/{{date}}/assembly/{sample}/{assembler}/{sample}.contigs.fasta".format(
            assembler=return_assembler(opt_sample), sample=opt_sample
        )
    )


def get_expanded_contigs(wildcards):
    return [
        get_contigs(wildcards, sample)
        for sample in get_samples_for_date(wildcards.date)
    ]


def get_read_counts(wildcards):
    return (
        "results/{date}/assembly/{assembler}/{sample}.log".format(
            assembler=return_assembler(wildcards.sample), **wildcards
        ),
    )


def get_bwa_index(wildcards):
    if wildcards.reference == "human" or wildcards.reference == "main+human":
        return rules.bwa_large_index.output
    else:
        return rules.bwa_index.output


def get_target_events(wildcards):
    if wildcards.reference == "main" or wildcards.clonality != "clonal":
        # calling variants against the wuhan reference or we are explicitly interested in subclonal as well
        return "SUBCLONAL_MINOR SUBCLONAL_MAJOR SUBCLONAL_HIGH CLONAL"
    else:
        # only keep clonal variants
        return "CLONAL"


def get_control_fdr_input(wildcards):
    if wildcards.reference == "main":
        return "results/{date}/filtered-calls/ref~{reference}/{sample}.{filter}.bcf"
    else:
        # If reference is not main, we are polishing an assembly.
        # Here, there is no need to structural variants or annotation based filtering.
        # Hence we directly take the output of varlociraptor call on the small variants.
        return get_pattern_by_technology(
            wildcards,
            illumina_pattern="results/{date}/calls/ref~{reference}/{sample}.small.bcf",
            ont_pattern="results/{date}/calls/ref~{reference}/{sample}.bcf",
            ion_torrent_pattern="results/{date}/calls/ref~{reference}/{sample}.small.bcf",
        )


def get_vembrane_expression(wildcards):
    return config["variant-calling"]["filters"][wildcards.filter]


def zip_expand(expand_string, zip_wildcard_1, zip_wildcard_2, expand_wildcard):
    """
    Zip by two wildcards and the expand the zip over another wildcard.
    expand_string must contain {zip1}, {zip2} and {{exp}}.
    """

    return sum(
        [
            expand(ele, exp=expand_wildcard)
            for ele in expand(
                expand_string,
                zip,
                zip1=zip_wildcard_1,
                zip2=zip_wildcard_2,
            )
        ],
        [],
    )


def get_quast_fastas(wildcards):
    if wildcards.stage == "unpolished":
        return "results/{date}/contigs/checked/{sample}.fasta"
    elif wildcards.stage == "polished":
        return "results/{date}/contigs/polished/{sample}.fasta"
    elif wildcards.stage == "masked/polished":
        return "results/{date}/contigs/masked/polished/{sample}.fasta"
    elif wildcards.stage == "pseudoassembly":
        return "results/{date}/contigs/pseudoassembled/{sample}.fasta"
    elif wildcards.stage == "masked/consensus":
        return "results/{date}/contigs/masked/consensus/{sample}.fasta"


def get_random_strain():
    with checkpoints.extract_strain_genomes_from_gisaid.get(
        date=BENCHMARK_DATE_WILDCARD
    ).output[0].open() as f:
        lines = f.read().splitlines()
        rnd_strain_path = random.choice(lines)
        strain = rnd_strain_path.replace(".fasta", "").split("/")[-1]
        return strain


def generate_mixtures(wildcards):
    if not config["mixtures"]["use_predefined_mixtures"]:
        no_mixtures = config["mixtures"]["no_mixtures"]
        no_strains = config["mixtures"]["no_strains"]
        mixture_list = []

        for mix in range(no_mixtures):

            fractions = [random.randint(1, 100) for _ in range(no_strains)]
            s = sum(fractions)
            fractions = [round(i / s * 100) for i in fractions]

            s = sum(fractions)
            if s != 100:
                fractions[-1] += 100 - s

            mixture = ""
            for frac in fractions:
                strain = get_random_strain()
                mixture += f"{MIXTURE_PART_INDICATOR}{strain}{MIXTURE_PERCENTAGE_INDICATOR}{frac}"

            mixture_list.append(mixture.replace(".", "-"))
    else:
        mixture_list = config["mixtures"]["predefined_mixtures"]
    return mixture_list


def get_mixture_results(wildcards):
    mixture_list = []

    with checkpoints.generate_mixtures.get().output[0].open() as f:
        for mix in f.read().splitlines():
            mixture_list.append(mix)

    if wildcards.caller == "pangolin":
        return expand(
            "results/benchmarking/tables/strain-calls/{prefix}{mixtures}.polished.strains.{caller}.csv",
            prefix=MIXTURE_PREFIX,
            caller=wildcards.caller,
            mixtures=mixture_list,
        )
    else:
        return expand(
            "results/benchmarking/tables/strain-calls/{prefix}{mixtures}.strains.{caller}.tsv",
            prefix=MIXTURE_PREFIX,
            caller=wildcards.caller,
            mixtures=mixture_list,
        )


def get_genome_fasta(wildcards):
    # mixtures sample, use provided (GISAID) genomes to generate mixtures of these
    if (
        MIXTURE_PART_INDICATOR in wildcards.accession
        and MIXTURE_PERCENTAGE_INDICATOR in wildcards.accession
    ):
        with checkpoints.extract_strain_genomes_from_gisaid.get(
            date=BENCHMARK_DATE_WILDCARD
        ).output[0].open() as f:
            acc, _ = wildcards.accession.split(MIXTURE_PERCENTAGE_INDICATOR)
            acc = acc.replace("-", ".").replace(MIXTURE_PART_INDICATOR, "")
            return "resources/genomes/{accession}.fasta".format(accession=acc)
    # read test sample
    if (
        READ_NUMBER_INDICATOR in wildcards.accession
        and READ_LENGTH_INDICATOR in wildcards.accession
    ):
        with checkpoints.extract_strain_genomes_from_gisaid.get(
            date=BENCHMARK_DATE_WILDCARD
        ).output[0].open() as f:
            acc, _ = wildcards.accession.split(READ_NUMBER_INDICATOR)
            acc = acc.replace(READ_TEST_PREFIX, "").replace("-", ".")
            return "resources/genomes/{accession}.fasta".format(accession=acc)
    # normal genome, download via entrez
    else:
        return "resources/genomes/{accession}.fasta".format(
            accession=wildcards.accession
        )


def no_reads(wildcards):
    max_reads = config["mixtures"]["max_reads"]
    if MIXTURE_PART_INDICATOR in wildcards.accession:
        _, fraction = wildcards.accession.split(MIXTURE_PERCENTAGE_INDICATOR)
        return round(int(fraction) * max_reads / 100)
    if READ_NUMBER_INDICATOR in wildcards.accession:
        _, no_reads = wildcards.accession.split(READ_NUMBER_INDICATOR)
        no_reads, _ = no_reads.split(READ_LENGTH_INDICATOR)
        return no_reads
    else:
        return max_reads


def length_read(wildcards):
    if READ_LENGTH_INDICATOR in wildcards.accession:
        _, length_state = wildcards.accession.split(READ_LENGTH_INDICATOR)
        length, _ = length_state.split(READ_STATE_INDICATOR)
        return length
    return 100


def get_strain(path_to_pangolin_call):
    pangolin_results = pd.read_csv(path_to_pangolin_call)
    return pangolin_results.loc[0]["lineage"]


def is_benchmark_data(sample):
    if (
        sample.startswith(BENCHMARK_PREFIX)
        or sample.startswith(NON_COV2_TEST_PREFIX)
        or sample.startswith(MIXTURE_PREFIX)
        or sample.startswith(READ_TEST_PREFIX)
    ):
        return True
    return False


def is_amplicon_data(sample):
    if is_benchmark_data(sample):
        # benchmark data, not amplicon based
        return False
    sample = pep.sample_table.loc[sample]
    try:
        return bool(int(sample["is_amplicon_data"]))
    except KeyError:
        return False


def any_sample_is_amplicon(wildcards):
    """Returns if any sample of the date is based on amplicon data."""
    return any(
        is_amplicon_data(sample) for sample in get_samples_for_date(wildcards.date)
    )


def get_list_of_amplicon_states(wildcards):
    return [True if is_amplicon_data(s) else False for s in get_samples()]


def get_list_of_amplicon_states_assembler(samples):
    return [return_assembler(s) for s in samples]


def get_varlociraptor_bias_flags(wildcards):
    if is_amplicon_data(wildcards.sample):
        # no bias detection possible
        return (
            "--omit-strand-bias --omit-read-orientation-bias --omit-read-position-bias"
        )
    return ""


def get_recal_input(wildcards):
    if is_amplicon_data(wildcards.sample):
        # do not mark duplicates
        return "results/{date}/mapped/ref~{reference}/{sample}.bam"
    # use BAM with marked duplicates
    return "results/{date}/dedup/ref~{reference}/{sample}.bam"


def get_depth_input(wildcards):
    # use clipped reads
    if is_amplicon_data(wildcards.sample):
        return get_pattern_by_technology(
            wildcards,
            illumina_pattern="results/{date}/read-sorted/pe~position/{sample}.hardclipped.bam",
            ont_pattern=expand(
                "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam",
                ref=config["preprocessing"]["amplicon-reference"],
            ),
            ion_torrent_pattern="results/{date}/read-sorted/se~position/{sample}.hardclipped.bam",
        )

    # use trimmed reads
    return "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam".format(
        ref=config["preprocessing"]["amplicon-reference"]
    )


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
        return config["preprocessing"]["kit-adapters"]


def get_final_assemblies(wildcards):
    if wildcards.assembly_type == "masked-assembly":
        return expand(
            "results/{{date}}/contigs/masked/polished/{sample}.fasta",
            sample=get_samples_for_date(wildcards.date),
        )
    elif wildcards.assembly_type == "pseudo-assembly":
        return get_list_of_expanded_patters_by_technology(
            wildcards,
            illumina_pattern="results/{{date}}/contigs/pseudoassembled/{sample}.fasta",
            ion_torrent_pattern="results/{{date}}/contigs/pseudoassembled/{sample}.fasta",
        )
    elif wildcards.assembly_type == "consensus-assembly":
        return get_list_of_expanded_patters_by_technology(
            wildcards,
            ont_pattern="results/{{date}}/contigs/masked/consensus/{sample}.fasta",
        )


def get_final_assemblies_identity(wildcards):
    if wildcards.assembly_type == "masked-assembly":
        return expand(
            "results/{{date}}/quast/masked/polished/{sample}/report.tsv",
            sample=get_samples_for_date(wildcards.date),
        )
    elif wildcards.assembly_type == "pseudo-assembly":
        return get_list_of_expanded_patters_by_technology(
            wildcards,
            illumina_pattern="results/{{date}}/quast/pseudoassembly/{sample}/report.tsv",
            ion_torrent_pattern="results/{{date}}/quast/pseudoassembly/{sample}/report.tsv",
        )
    elif wildcards.assembly_type == "consensus-assembly":
        return get_list_of_expanded_patters_by_technology(
            wildcards,
            ont_pattern="results/{{date}}/quast/masked/consensus/{sample}/report.tsv",
        )


def get_checkpoints_for_overview_table(wildcards):
    assembly_types = ["masked-assembly"]

    all_samples_for_date = get_samples_for_date(wildcards.date)

    if any(has_pseudo_assembly(None, sample) for sample in all_samples_for_date):
        assembly_types.append("pseudo-assembly")

    if any(has_consensus_assembly(None, sample) for sample in all_samples_for_date):
        assembly_types.append("consensus-assembly")

    return expand(
        "results/{{date}}/tables/quality-filter/{assembly_type}.txt",
        assembly_type=assembly_types,
    )


def load_filtered_samples(wildcards, assembly_type):
    with checkpoints.quality_filter.get(
        date=wildcards.date, assembly_type=assembly_type
    ).output.passed_filter.open() as f:
        try:
            return pd.read_csv(f, squeeze=True, header=None).astype(str).to_list()
        except pd.errors.EmptyDataError:
            return []


def get_assemblies_for_submission(wildcards, agg_type):
    if "sample" in wildcards:
        if wildcards.sample.startswith(READ_TEST_PREFIX):
            _, state = wildcards.sample.split(READ_STATE_INDICATOR)
            if state == "contig":
                return "results/{date}/tables/largest_contig/{sample}.fasta"
            elif state == "scaffold":
                return "results/{date}/contigs/ordered/{sample}.fasta"
            elif state == "polished_scaffold":
                return "results/{date}/contigs/polished/{sample}.fasta"
            elif state == "pseudo":
                return "results/{date}/contigs/pseudoassembled/{sample}.fasta"

    if wildcards.date != BENCHMARK_DATE_WILDCARD:

        all_samples_for_date = get_samples_for_date(wildcards.date)

        masked_samples = load_filtered_samples(wildcards, "masked-assembly")
        pseudo_samples = (
            load_filtered_samples(wildcards, "pseudo-assembly")
            if any(has_pseudo_assembly(None, sample) for sample in all_samples_for_date)
            else []
        )
        consensus_samples = (
            load_filtered_samples(wildcards, "consensus-assembly")
            if any(
                has_consensus_assembly(None, sample) for sample in all_samples_for_date
            )
            else []
        )

    # for testing of pangolin don't create pseudo-assembly
    else:
        masked_samples = [wildcards.sample]

    normal_assembly_pattern = "results/{{date}}/contigs/masked/polished/{sample}.fasta"
    pseudo_assembly_pattern = "results/{{date}}/contigs/pseudoassembled/{sample}.fasta"
    consensus_assembly_pattern = (
        "results/{{date}}/contigs/masked/consensus/{sample}.fasta"
    )

    # get accepted samples for rki submission
    if agg_type == "accepted samples" or agg_type == "accepted samples technology":
        selected_assemblies = []
        technolgy = []
        unqiue_samples = set()

        if len(set(masked_samples)) > 0:
            unqiue_samples.update(masked_samples)
        if len(set(pseudo_samples)) > 0:
            unqiue_samples.update(pseudo_samples)
        if len(set(consensus_samples)) > 0:
            unqiue_samples.update(consensus_samples)

        if len(unqiue_samples) == 0:
            # No sequences are passing the quality filter for date.
            # Return dummy path as an indicator
            return "resources/genomes/main.fasta"

        for sample in unqiue_samples:
            if sample in masked_samples:
                selected_assemblies.append(
                    normal_assembly_pattern.format(sample=sample)
                )
            elif sample in pseudo_samples:
                selected_assemblies.append(
                    pseudo_assembly_pattern.format(sample=sample)
                )
            elif sample in consensus_samples:
                selected_assemblies.append(
                    consensus_assembly_pattern.format(sample=sample)
                )

            technolgy.append(
                get_pattern_by_technology(
                    wildcards,
                    sample=sample,
                    illumina_pattern="ILLUMINA",
                    ont_pattern="OXFORD_NANOPORE",
                    ion_torrent_pattern="ION_TORRENT",
                )
            )

        if agg_type == "accepted samples":
            return selected_assemblies
        elif agg_type == "accepted samples technology":
            return technolgy

    # for the pangolin call
    elif agg_type == "single sample":
        if wildcards.sample in masked_samples:
            return "results/{date}/contigs/polished/{sample}.fasta"
        elif wildcards.sample in pseudo_samples:
            return "results/{date}/contigs/pseudoassembled/{sample}.fasta"
        elif wildcards.sample in consensus_samples:
            return "results/{date}/contigs/consensus/{sample}.fasta"
        # for not accepted samples call on the polished-contigs
        else:
            return "results/{date}/contigs/polished/{sample}.fasta"

    # for the qc report
    elif agg_type == "all samples":
        assembly_type_used = []
        for sample in all_samples_for_date:
            if sample in masked_samples:
                assembly_type_used.append(f"{sample},normal")
            elif sample in pseudo_samples:
                assembly_type_used.append(f"{sample},pseudo")
            elif sample in consensus_samples:
                assembly_type_used.append(f"{sample},consensus")
            else:
                assembly_type_used.append(f"{sample},not-accepted")
        return assembly_type_used

    return inner


def get_output_dir(wildcards, output):
    return os.path.dirname(output[0])


def expand_samples_by_func(paths, func, **kwargs):
    def inner(wildcards):
        return expand(
            paths,
            sample=func(wildcards.date),
            **kwargs,
        )

    return inner


def expand_samples_for_date(paths, **kwargs):
    return expand_samples_by_func(paths, get_samples_for_date, **kwargs)


def get_input_plotting_primer_clipping(wildcards, stage, suffix=""):
    """Returns list of unclipped bam files for a date. Used for visualizing the primer clipping."""
    return get_list_of_expanded_patters_by_technology(
        wildcards,
        return_only_amplicon_samples=True,
        illumina_pattern=f"results/{{{{date}}}}/read-sorted/pe~position/{{sample}}.{stage}.bam{suffix}",
        ont_pattern=f"results/{{{{date}}}}/read-sorted/se~position/{{sample}}.{stage}.bam{suffix}",
        ion_torrent_pattern=f"results/{{{{date}}}}/read-sorted/se~position/{{sample}}.{stage}.bam{suffix}",
    )


def get_fallbacks_for_report(fallback_type):
    """Returns path to the fallback sequences. The "main.fasta" is returned as an indicator that no fallback sequences is created."""

    def inner(wildcards):
        samples = get_samples_for_date(wildcards.date)

        if fallback_type == "pseudo":
            path = "results/{{date}}/contigs/pseudoassembled/{sample}.fasta"
            return [
                path.format(sample=sample)
                if has_pseudo_assembly(None, sample)
                else "resources/genomes/main.fasta"
                for sample in get_samples_for_date(wildcards.date)
            ]

        elif fallback_type == "consensus":
            path = "results/{{date}}/contigs/masked/consensus/{sample}.fasta"
            return [
                path.format(sample=sample)
                if has_consensus_assembly(None, sample)
                else "resources/genomes/main.fasta"
                for sample in samples
            ]

        raise NotImplementedError(f'No fallback for "{fallback_type}" found.')

    return inner


def get_pattern_by_technology(
    wildcards,
    illumina_pattern=None,
    ont_pattern=None,
    ion_torrent_pattern=None,
    sample=None,
):
    """Returns the given pattern, depending on the sequencing technology used for the sample."""
    if sample is None:
        if is_illumina(wildcards):
            return illumina_pattern
        elif is_ont(wildcards):
            return ont_pattern
        elif is_ion_torrent(wildcards):
            return ion_torrent_pattern

    if is_illumina(None, sample):
        return illumina_pattern
    elif is_ont(None, sample):
        return ont_pattern
    elif is_ion_torrent(None, sample):
        return ion_torrent_pattern

    raise NotImplementedError(
        f'The technolgy listed for sample "{wildcards.sample}" is not supported.'
    )


def format_patterns(input_patterns, sample, formated_patterns):
    """Add the sample to the given pattern, depending if its a string or a list."""

    if isinstance(input_patterns, str):
        formated_patterns.append(input_patterns.format(sample=sample))
    elif isinstance(input_patterns, list):
        [
            formated_patterns.append(pattern.format(sample=sample))
            for pattern in input_patterns
        ]
    else:
        raise TypeError()

    return formated_patterns


def get_list_of_expanded_patters_by_technology(
    wildcards,
    illumina_pattern=None,
    ont_pattern=None,
    ion_torrent_pattern=None,
    return_only_amplicon_samples=False,
):
    """Returns an aggregate list of given patterns, depending on their sequences technology. Formats the {sample} wildcard."""
    patterns = []

    samples = get_samples_for_date(wildcards.date)

    if return_only_amplicon_samples:
        for sample in samples:
            if (
                illumina_pattern is not None
                and is_illumina(None, sample)
                and is_amplicon_data(sample)
            ):
                patterns = format_patterns(illumina_pattern, sample, patterns)
            elif (
                ont_pattern is not None
                and is_ont(None, sample)
                and is_amplicon_data(sample)
            ):
                patterns = format_patterns(ont_pattern, sample, patterns)
            elif (
                ion_torrent_pattern is not None
                and is_ion_torrent(None, sample)
                and is_amplicon_data(sample)
            ):
                patterns = format_patterns(ion_torrent_pattern, sample, patterns)
        return patterns

    for sample in samples:
        if illumina_pattern is not None and is_illumina(None, sample):
            patterns = format_patterns(illumina_pattern, sample, patterns)
        elif ont_pattern is not None and is_ont(None, sample):
            patterns = format_patterns(ont_pattern, sample, patterns)
        elif ion_torrent_pattern is not None and is_ion_torrent(None, sample):
            patterns = format_patterns(ion_torrent_pattern, sample, patterns)
    return patterns


def get_raw_reads_counts(wildcards):
    """Returns paths of files to be parsed by the overview table rule for the raw reads counts."""
    return get_list_of_expanded_patters_by_technology(
        wildcards,
        illumina_pattern="results/{{date}}/trimmed/fastp-pe/{sample}.fastp.json",
        ont_pattern="results/{{date}}/tables/fastq-read-counts/raw~{sample}.txt",
        ion_torrent_pattern="results/{{date}}/trimmed/fastp-se/{sample}.fastp.json",
    )


def get_trimmed_reads_counts(wildcards):
    """Return paths of files to be parsed by the overview table rule for the trimmed reads counts."""
    return get_list_of_expanded_patters_by_technology(
        wildcards,
        illumina_pattern="results/{{date}}/trimmed/fastp-pe/{sample}.fastp.json",
        ont_pattern="results/{{date}}/tables/fastq-read-counts/trimmed~{sample}.txt",
        ion_torrent_pattern="results/{{date}}/trimmed/fastp-se/{sample}.fastp.json",
    )


def get_fastp_results(wildcards):
    """Returns paths of files to aggregate the fastp results for the multiqc rule."""
    # fastp is only used on Illumina and Ion Torrent data
    return get_list_of_expanded_patters_by_technology(
        wildcards,
        illumina_pattern="results/{{date}}/trimmed/fastp-pe/{sample}.fastp.json",
        ion_torrent_pattern="results/{{date}}/trimmed/fastp-se/{sample}.fastp.json",
    )


def get_vep_args(wildcards, input):
    return (
        "--vcf_info_field ANN --hgvsg --hgvs --synonyms {synonyms} "
        "--custom {input.problematic},,vcf,exact,0,"
    ).format(input=input, synonyms=get_resource("synonyms.txt"))


def get_samples_for_assembler_comparison(paths):
    return zip_expand(
        expand_string=paths,
        zip_wildcard_1=get_illumina_samples_dates(),
        zip_wildcard_2=get_illumina_samples(),
        expand_wildcard=config["assemblers_for_comparison"],
    )


def get_illumina_samples():
    samples = pep.sample_table
    return samples.loc[samples["technology"] == ILLUMINA]["sample_name"].values


def get_illumina_samples_dates():
    samples = pep.sample_table
    return samples.loc[samples["technology"] == ILLUMINA]["date"].values


def get_megahit_preset(wildcards):
    if wildcards.preset == "std":
        return ""
    else:
        return f"--preset {wildcards.preset}"


def get_lineage_by_accession(wildcards):
    return list(config["strain-calling"]["lineage-references"].keys())[
        list(config["strain-calling"]["lineage-references"].values()).index(
            wildcards.accession
        )
    ]


def get_include_flag(sample):
    try:
        samples = pep.sample_table
        samples.dropna(subset=["include_in_high_genome_summary"], inplace=True)
        return samples.loc[sample]["include_in_high_genome_summary"]
    # if there is no include_in_high_genome_summary in the
    # samples.csvdefined, always include the sample
    except KeyError:
        return 1


def get_include_flag_for_date(wildcards):
    return [
        get_include_flag(sample)
        for sample in get_assemblies_for_submission(wildcards, "accepted samples")
    ]


def get_artic_primer(wildcards):
    # TODO add more _adapters.py (not preferred) or
    # add a script to generate them from a link to a bed file.
    # The bed file can be found in the artic repo. Related to #356
    return "resources/ARTIC_v{}_adapters.py".format(
        config["preprocessing"]["artic-primer-version"]
    )


def get_trimmed_reads(wildcards):
    """Returns paths of files of the trimmed reads for parsing by kraken."""
    return get_list_of_expanded_patters_by_technology(
        wildcards,
        illumina_pattern=expand(
            "results/{{{{date}}}}/trimmed/fastp-pe/{{sample}}.{read}.fastq.gz",
            read=[1, 2],
        ),
        ont_pattern="results/{{date}}/trimmed/porechop/adapter_barcode_trimming/{sample}.fastq.gz",
        ion_torrent_pattern="results/{{date}}/trimmed/fastp-se/{sample}.fastq.gz",
    )


def get_kraken_output(wildcards):
    """Returns the output of kraken on the raw reads, depend on sequencing technology."""
    return get_list_of_expanded_patters_by_technology(
        wildcards,
        illumina_pattern="results/{date}/species-diversity/pe/{{sample}}/{{sample}}.uncleaned.kreport2".format(
            **wildcards
        ),
        ont_pattern="results/{date}/species-diversity/se/{{sample}}/{{sample}}.uncleaned.kreport2".format(
            **wildcards
        ),
        ion_torrent_pattern="results/{date}/species-diversity/se/{{sample}}/{{sample}}.uncleaned.kreport2".format(
            **wildcards
        ),
    )


def get_kraken_output_after_filtering(wildcards):
    """Returns the output of kraken on the filtered reads, depend on sequencing technology."""
    return get_list_of_expanded_patters_by_technology(
        wildcards,
        illumina_pattern="results/{date}/species-diversity-nonhuman/pe/{{sample}}/{{sample}}.cleaned.kreport2".format(
            **wildcards
        ),
        ont_pattern="results/{date}/species-diversity-nonhuman/se/{{sample}}/{{sample}}.cleaned.kreport2".format(
            **wildcards
        ),
        ion_torrent_pattern="results/{date}/species-diversity-nonhuman/se/{{sample}}/{{sample}}.cleaned.kreport2".format(
            **wildcards
        ),
    )


def get_read_calls(wildcards):
    with checkpoints.select_random_lineages.get(date=BENCHMARK_DATE_WILDCARD).output[
        0
    ].open() as f:
        lineages = f.read().splitlines()

    return expand(
        "results/benchmarking/tables/collected_lineage_calls_on_{lineage}_{number}_{length}.tsv",
        lineage=lineages,
        number=config["read_lineage_call"]["number_of_reads"],
        length=config["read_lineage_call"]["length_of_reads"],
    )


def get_first_line(path):
    with open(path) as f:
        return f.readline().strip()


def get_kallisto_quant_extra(wildcards, input):
    if is_for_testing():
        return get_if_testing("--single --fragment-length 250 --sd 47301")

    return (
        f"--single --fragment-length {get_first_line(input.fragment_length)} --sd {get_first_line(input.standard_deviation)}"
        if is_single_end(wildcards)
        else "",
    )


def get_kallisto_quant_input(wildcards):
    if is_single_end(wildcards):
        return {
            "fastq": get_reads_after_qc(wildcards),
            "index": "results/{date}/kallisto/strain-genomes.idx",
            "fragment_length": "results/{date}/tables/avg_read_length/{sample}.txt",
            "standard_deviation": "results/{date}/tables/standard_deviation/{sample}.txt",
        }
    return {
        "fastq": get_reads_after_qc(wildcards),
        "index": "results/{date}/kallisto/strain-genomes.idx",
    }


def is_for_testing():
    return bool(config.get("testing", {}))


def get_if_testing(string):
    return string if is_for_testing() else ""


def get_reads_by_stage(wildcards):
    if wildcards.stage == "raw":
        return get_fastqs(wildcards)
    elif wildcards.stage == "trimmed":
        return "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fastq"
    elif wildcards.stage == "clipped":
        return "results/{date}/trimmed/porechop/primer_clipped/{sample}.fastq"
    elif wildcards.stage == "filtered":
        return "results/{date}/trimmed/nanofilt/{sample}.fastq"


def get_polished_sequence(wildcards):
    """Returns path to polished sequences, depend on sequencing technology."""
    return get_pattern_by_technology(
        wildcards,
        illumina_pattern="results/{date}/polishing/bcftools-illumina/{sample}.fasta",
        ont_pattern="results/{date}/polishing/medaka/{sample}/{sample}.fasta",
        ion_torrent_pattern="results/{date}/polishing/bcftools-illumina/{sample}.fasta",
    )


def get_fallback_sequence(wildcards):
    """Returns path to fallback sequences, depend on sequencing technology."""
    return get_pattern_by_technology(
        wildcards,
        illumina_pattern="results/{date}/contigs/pseudoassembled/{sample}.fasta",
        ont_pattern="results/{date}/contigs/consensus/{sample}.fasta",
        ion_torrent_pattern="results/{date}/contigs/pseudoassembled/{sample}.fasta",
    )


def get_varrange(wildcards):
    return get_pattern_by_technology(
        wildcards,
        illumina_pattern=ILLUMINA_VARRANGE,
        ont_pattern=ONT_VARRANGE,
        ion_torrent_pattern=ION_VARRANGE,
    )


def get_if_any_consensus_assembly(path):
    """Returns the samples for which consensus-assemblies should be created."""

    def inner(wildcards):
        if any(
            has_consensus_assembly(None, sample)
            for sample in get_samples_for_date(wildcards.date)
        ):
            return path
        return "resources/genomes/main.fasta"

    return inner


def get_if_any_pseudo_assembly(path):
    """Returns the samples for which pseudo-assemblies should be created."""

    def inner(wildcards):
        if any(
            has_pseudo_assembly(None, sample)
            for sample in get_samples_for_date(wildcards.date)
        ):
            return path
        return "resources/genomes/main.fasta"

    return inner


def get_seq_type(wildcards):
    """Returns the sequencing type used for the samples of a given date."""
    # see: https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Anleitung-Bereitstellung-Sequenzdaten.pdf?__blob=publicationFile
    return get_list_of_expanded_patters_by_technology(
        wildcards,
        illumina_pattern="ILLUMINA",
        ont_pattern="OXFORD_NANOPORE",
        ion_torrent_pattern="ION_TORRENT",
    )


def get_samtools_sort_input(wildcards):
    """Returns the input for the samtools sort rule."""
    if wildcards.stage == "initial":
        return expand(
            "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam",
            ref=config["preprocessing"]["amplicon-reference"],
        )
    elif wildcards.stage == "hardclipped":
        return (
            "results/{date}/read-clipping/hardclipped/{read_type}/{sample}/{sample}.bam"
        )

    raise NotImplementedError(f"Sorting for {wildcards.stage} not supported.")


def get_pangolin_input(wildcards):
    if wildcards.stage == "scaffold":
        return "results/{date}/contigs/ordered/{sample}.fasta"
    elif wildcards.stage == "polished":
        return "results/{date}/contigs/polished/{sample}.fasta"
    elif wildcards.stage == "masked-polished":
        return "results/{date}/contigs/masked/polished/{sample}.fasta"
    elif wildcards.stage == "pseudo":
        return "results/{date}/contigs/pseudoassembled/{sample}.fasta"
    elif wildcards.stage == "consensus":
        return "results/{date}/consensus/bcftools/{sample}.fasta"
    elif wildcards.stage == "masked-consensus":
        return "results/{date}/contigs/masked/consensus/{sample}.fasta"


def get_pangolin_stage_by_technolgy(sample):
    if has_pseudo_assembly(None, sample):
        return ["scaffold", "polished", "masked-polished", "pseudo"]
    elif has_consensus_assembly(None, sample):
        return [
            "scaffold",
            "polished",
            "masked-polished",
            "consensus",
            "masked-consensus",
        ]

    raise NotImplementedError(f"No pangolin stages for technology {technology} found.")


def get_aggregated_pangolin_calls(wildcards, return_list="paths"):
    samples = get_samples_for_date(wildcards.date)

    pangolin_pattern = (
        "results/{date}/tables/strain-calls/{sample}.{stage}.strains.pangolin.csv"
    )
    expanded_patterns = []

    for sample in samples:

        stage_wildcards = get_pattern_by_technology(
            wildcards,
            sample=sample,
            illumina_pattern=get_pangolin_stage_by_technolgy(sample),
            ont_pattern=get_pangolin_stage_by_technolgy(sample),
            ion_torrent_pattern=get_pangolin_stage_by_technolgy(sample),
        )

        for stage in stage_wildcards:
            if return_list == "paths":
                expanded_patterns.append(
                    pangolin_pattern.format(
                        stage=stage, sample=sample, date=wildcards.date
                    )
                )
            elif return_list == "stages":
                expanded_patterns.append(stage)
            elif return_list == "samples":
                expanded_patterns.append(sample)
            else:
                raise NameError(f"return_list {return_list} not recognized.")

    return expanded_patterns


def get_pangolin_for_report(wildcards):
    paths = []

    path = "results/{date}/tables/strain-calls/{sample}.{stage}.strains.pangolin.csv"

    for entry in get_assemblies_for_submission(wildcards, "all samples"):
        sample, assembly = entry.split(",")
        if assembly == "normal":
            pango_stage = "polished"
        elif assembly == "pseudo":
            pango_stage = "pseudo"
        elif assembly == "consensus":
            pango_stage = "consensus"
        elif assembly == "not-accepted":
            pango_stage = "polished"

        paths.append(path.format(sample=sample, date=wildcards.date, stage=pango_stage))

    return paths


wildcard_constraints:
    sample="[^/.]+",
    vartype="|".join(VARTYPES),
    clonality="subclonal|clonal",
    filter="|".join(
        list(map(re.escape, config["variant-calling"]["filters"])) + ["nofilter"]
    ),
    varrange="structural|small|homopolymer-medaka|homopolymer-longshot",
