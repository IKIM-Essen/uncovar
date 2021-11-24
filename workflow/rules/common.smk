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

ILLUMINA = "illumina"
ONT = "ont"

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


def is_sample_illumina(sample):
    return get_technology_by_sample(sample) == ILLUMINA


def is_sample_ont(sample):
    return get_technology_by_sample(sample) == ONT


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
    if is_ont(wildcards):
        return pep.sample_table.loc[wildcards.sample][["fq1"]]
    elif is_illumina(wildcards):
        return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]


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
        "results/benchmarking/tables/strain-calls/non-cov2-{accession}.strains.pangolin.csv"
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
        elif wildcards.reference == config["adapters"]["amplicon-reference"]:
            # return reference genome of amplicon primers
            return "resources/genomes/{reference}.fasta{suffix}".format(
                reference=config["adapters"]["amplicon-reference"], suffix=suffix
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
        if is_illumina(wildcards):
            return expand(
                "results/{date}/trimmed/{sample}.{read}.fastq.gz",
                read=[1, 2],
                **wildcards,
            )
        elif is_ont(wildcards):
            return expand(
                "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
                **wildcards,
            )

    # theses reads are used to generate the bam file for the BAMclipper and the coverage plot of the main reference
    elif wildcards.reference == config["adapters"]["amplicon-reference"]:
        return get_non_human_reads(wildcards)

    # aligments to other references (e.g. the covid reference genome),
    # are done with reads, which have undergone the quality control process
    else:
        return get_reads_after_qc(wildcards)


def get_non_human_reads(wildcards):
    if is_illumina(wildcards):
        return expand(
            "results/{date}/nonhuman-reads/pe/{sample}.{read}.fastq.gz",
            date=wildcards.date,
            read=[1, 2],
            sample=wildcards.sample,
        )
    elif is_ont(wildcards):
        return expand(
            "results/{date}/nonhuman-reads/se/{sample}.fastq.gz",
            date=wildcards.date,
            sample=wildcards.sample,
        )


def get_reads_after_qc(wildcards, read="both"):
    if is_amplicon_data(wildcards.sample) and is_ont(wildcards):
        pattern = [
            "results/{date}/nonhuman-reads/se/{sample}.fastq.gz".format(**wildcards)
        ]
    elif not is_amplicon_data(wildcards.sample) and is_ont(wildcards):
        raise NotImplementedError(
            "UnCoVer currently does not support non amplicon based ONT data"
        )
    elif is_amplicon_data(wildcards.sample) and is_illumina(wildcards):
        pattern = expand(
            "results/{date}/clipped-reads/{sample}.{read}.fastq.gz",
            date=wildcards.date,
            read=[1, 2],
            sample=wildcards.sample,
        )
    elif not is_amplicon_data(wildcards.sample) and is_illumina(wildcards):
        pattern = expand(
            "results/{date}/nonhuman-reads/pe/{sample}.{read}.fastq.gz",
            date=wildcards.date,
            read=[1, 2],
            sample=wildcards.sample,
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
    if is_amplicon_data(sample) and is_illumina(None, sample):
        return config["assembly"]["amplicon"]
    elif not is_amplicon_data(sample) and is_illumina(None, sample):
        return config["assembly"]["shotgun"]
    if is_amplicon_data(sample) and is_ont(None, sample):
        raise NotImplementedError(
            f"No amplicon ONT assembler option found for sample {sample}"
        )
    elif not is_amplicon_data(sample) and is_ont(None, sample):
        raise NotImplementedError(
            f"No non-amplicon ONT assembler option found for sample {sample}"
        )


def get_contigs(wildcards, opt_sample=None):

    if "sample" in wildcards.keys():
        if is_ont(wildcards):
            return "results/{date}/assembly/{sample}/spades_se/{sample}.contigs.fasta"

        elif is_illumina(wildcards):
            return "results/{{date}}/assembly/{{sample}}/{assembler}/{{sample}}.contigs.fasta".format(
                assembler=return_assembler(wildcards.sample)
            )

    # wildcards is only sample name
    if is_ont(None, opt_sample):
        return "results/{{date}}/assembly/{sample}/spades_se/{sample}.contigs.fasta".format(
            sample=opt_sample
        )

    elif is_illumina(None, opt_sample):
        return "results/{{date}}/assembly/{sample}/{assembler}/{sample}.contigs.fasta".format(
            assembler=return_assembler(opt_sample), sample=opt_sample
        )

    raise NotImplementedError("No assembler found.")


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


def get_filter_odds_input(wildcards):
    if wildcards.reference == "main":
        return "results/{date}/filtered-calls/ref~{reference}/{sample}.{filter}.bcf"
    else:
        # If reference is not main, we are polishing an assembly.
        # Here, there is no need to structural variants or annotation based filtering.
        # Hence we directly take the output of varlociraptor call on the small variants.
        return "results/{date}/calls/ref~{reference}/{sample}.small.bcf"


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
        return get_contigs(wildcards)
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
            "results/benchmarking/tables/strain-calls/{prefix}{mixtures}.strains.{caller}.csv",
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


def get_samples_for_date_for_illumina_amplicon(date):
    return [
        s
        for s in get_samples_for_date(date)
        if (is_amplicon_data(s) and is_illumina(None, s))
    ]


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
    if is_illumina(wildcards) and is_amplicon_data(wildcards.sample):
        return "results/{date}/clipped-reads/{sample}.primerclipped.bam"

    elif is_ont(wildcards) and is_amplicon_data(wildcards.sample):
        return expand(
            "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam",
            ref=config["adapters"]["amplicon-reference"],
        )

    # use trimmed reads
    return "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam".format(
        ref=config["adapters"]["amplicon-reference"]
    )


def get_adapters(wildcards):
    if is_illumina(wildcards) and is_amplicon_data(wildcards.sample):
        return config["adapters"]["illumina-amplicon"]
    elif is_illumina(wildcards) and not is_amplicon_data(wildcards.sample):
        return config["adapters"]["illumina-shotgun"]
    elif is_ont(wildcards) and is_amplicon_data(wildcards.sample):
        raise NotImplementedError(
            "No adapters implemented for amplicon data generated with ONT technology"
        )
    elif is_ont(wildcards) and not is_amplicon_data(wildcards.sample):
        raise NotImplementedError(
            "No adapters implemented for shotgun data generated with ONT technology"
        )


def get_final_assemblies(wildcards):
    all_samples = get_samples_for_date(wildcards.date)
    illumina_samples = [sample for sample in all_samples if is_illumina(None, sample) ]
    ont_samples = [sample for sample in all_samples if is_ont(None, sample) ]

    if wildcards.assembly_type == "masked-assembly":
        return expand("results/{{date}}/contigs/masked/polished/{sample}.fasta", sample=all_samples)
    elif wildcards.assembly_type == "pseudo-assembly":
        return expand("results/{{date}}/contigs/pseudoassembled/{sample}.fasta", sample=illumina_samples)
    elif wildcards.assembly_type == "consensus-assembly":
        return expand("results/{{date}}/contigs/masked/consensus/{sample}.fasta", sample=ont_samples)


def get_final_assemblies_identity(wildcards):
    all_samples = get_samples_for_date(wildcards.date)
    illumina_samples = [sample for sample in all_samples if is_illumina(None, sample) ]
    ont_samples = [sample for sample in all_samples if is_ont(None, sample) ]

    if wildcards.assembly_type == "masked-assembly":
        return expand("results/{{date}}/quast/masked/polished/{sample}/report.tsv", sample=all_samples)
    elif wildcards.assembly_type == "pseudo-assembly":
        return expand("results/{{date}}/quast/pseudoassembly/{sample}/report.tsv", sample=illumina_samples)
    elif wildcards.assembly_type == "consensus-assembly":
        return expand("results/{{date}}/quast/masked/consensus/{sample}/report.tsv", sample=ont_samples)


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
            if any(is_illumina(None, sample) for sample in all_samples_for_date)
            else []
        )
        consensus_samples = (
            load_filtered_samples(wildcards, "consensus-assembly")
            if any(is_ont(None, sample) for sample in all_samples_for_date)
            else []
        )

        print("masked_samples", masked_samples)
        print("pseudo_samples", pseudo_samples)
        print("consensus_samples", consensus_samples)

    # for testing of pangolin don't create pseudo-assembly
    else:
        masked_samples = [wildcards.sample]

    normal_assembly_pattern = "results/{{date}}/contigs/masked/polished/{sample}.fasta"
    pseudo_assembly_pattern = "results/{{date}}/contigs/pseudoassembled/{sample}.fasta"
    consensus_assembly_pattern = (
        "results/{{date}}/contigs/masked/consensus/{sample}.fasta"
    )

    # get accepted samples for rki submission
    if agg_type == "accepted samples":
        selected_assemblies = []
        unqiue_samples = set()

        if len(set(masked_samples)) > 0:
            unqiue_samples.update(masked_samples)
        if len(set(pseudo_samples)) > 0:
            unqiue_samples.update(pseudo_samples)
        if len(set(consensus_samples)) > 0:
            unqiue_samples.update(consensus_samples)

        if len(unqiue_samples) == 0:
            raise NotImplementedError("No sequences are passing the quality filter.")

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

        return selected_assemblies

    # for the pangolin call
    elif agg_type == "single sample":
        if wildcards.sample in masked_samples:
            return "results/{date}/contigs/masked/polished/{sample}.fasta"
        elif wildcards.sample in pseudo_samples:
            return "results/{date}/contigs/pseudoassembled/{sample}.fasta"
        elif wildcards.sample in consensus_samples:
            return "results/{date}/contigs/masked/consensus/{sample}.fasta"
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


def expand_samples_for_date_amplicon(paths, **kwargs):
    return expand_samples_by_func(
        paths, get_samples_for_date_for_illumina_amplicon, **kwargs
    )


def get_for_report_if_illumina_sample(path):
    def inner(wildcards):
        return [
            path.format(sample=sample)
            if is_illumina(None, sample)
            else "resources/genomes/main.fasta"
            for sample in get_samples_for_date(wildcards.date)
        ]

    return inner


def get_for_report_if_ont_sample(path):
    def inner(wildcards):
        return [
            path.format(sample=sample)
            if is_ont(None, sample)
            else "resources/genomes/main.fasta"
            for sample in get_samples_for_date(wildcards.date)
        ]

    return inner


def get_raw_reads_counts(wildcards):
    return [
        "results/{{date}}/trimmed/{sample}.fastp.json".format(sample=sample)
        if is_illumina(None, sample)
        else "results/{{date}}/tables/fastq-read-counts/raw~{sample}.txt".format(
            sample=sample
        )
        for sample in get_samples_for_date(wildcards.date)
    ]


def get_trimmed_reads_counts(wildcards):
    return [
        "results/{{date}}/trimmed/{sample}.fastp.json".format(sample=sample)
        if is_illumina(None, sample)
        else "results/{{date}}/tables/fastq-read-counts/trimmed~{sample}.txt".format(
            sample=sample
        )
        for sample in get_samples_for_date(wildcards.date)
    ]


def get_fastp_results(wildcards, **kwargs):
    # fastp is only used on illumina data
    return expand(
        "results/{{date}}/trimmed/{sample}.fastp.json",
        sample=[
            sample
            for sample in get_samples_for_date(wildcards.date)
            if is_illumina(None, sample)
        ],
    )


def get_vep_args(wildcards, input):
    return (
        "--vcf_info_field ANN --hgvsg --hgvs --synonyms {synonyms} "
        "--custom {input.problematic},,vcf,exact,0,"
    ).format(input=input, synonyms=get_resource("synonyms.txt"))


def get_samples_for_assembler_comparison(paths):
    return zip_expand(
        paths,
        get_illumina_samples_dates(),
        get_illumina_samples(),
        config["assemblers_for_comparison"],
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


def get_artic_primer(wildcards):
    # TODO: add more _adapters.py (not preferred) or
    # add a script to generate them from a link to a bed file.
    # The bed file can be found in the artic repo
    return "resources/ARTIC_v{}_adapters.py".format(
        config["adapters"]["artic-primer-version"]
    )


def get_trimmed_reads(wildcards):
    if is_illumina(wildcards):
        return expand(
            "results/{{date}}/trimmed/{{sample}}.{read}.fastq.gz", read=[1, 2]
        )
    elif is_ont(wildcards):
        return (
            "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fastq.gz"
        )


def get_kraken_output(wildcards):
    samples = get_samples_for_date(wildcards.date)

    illumina_pattern = (
        "results/{date}/species-diversity/pe/{sample}/{sample}.uncleaned.kreport2"
    )
    ont_pattern = (
        "results/{date}/species-diversity/se/{sample}/{sample}.uncleaned.kreport2"
    )

    return [
        illumina_pattern.format(sample=sample, **wildcards)
        if is_illumina(None, sample)
        else ont_pattern.format(sample=sample, **wildcards)
        for sample in samples
    ]


def get_kraken_output_after_filtering(wildcards):
    samples = get_samples_for_date(wildcards.date)

    illumina_pattern = "results/{date}/species-diversity-nonhuman/pe/{sample}/{sample}.cleaned.kreport2"
    ont_pattern = "results/{date}/species-diversity-nonhuman/se/{sample}/{sample}.cleaned.kreport2"

    return [
        illumina_pattern.format(sample=sample, **wildcards)
        if is_illumina(None, sample)
        else ont_pattern.format(sample=sample, **wildcards)
        for sample in samples
    ]


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
        if is_ont(wildcards)
        else "",
    )


def get_path_if_ont(paths):
    def inner(wildcards):
        if is_ont(wildcards):
            return paths
        return [""]

    return inner


def get_kallisto_quant_input(wildcards):
    if is_ont(wildcards):
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
    if is_illumina(wildcards):
        return "results/{date}/polishing/bcftools-illumina/{sample}.fasta"
    elif is_ont(wildcards):
        return "results/{date}/polishing/medaka/{sample}/{sample}.fasta"


def get_varrange(wildcards):
    if is_ont(wildcards):
        return ["homopolymer"]
    elif is_illumina(wildcards):
        return ["small", "structural"]


def get_if_any_sample_is_ont(path):
    def inner(wildcards):
        if any(is_ont(None, sample) for sample in get_samples_for_date(wildcards.date)):
            return path
        return "resources/genomes/main.fasta"

    return inner


def get_if_any_sample_is_illumina(path):
    def inner(wildcards):
        if any(
            is_illumina(None, sample) for sample in get_samples_for_date(wildcards.date)
        ):
            return path
        return "resources/genomes/main.fasta"

    return inner

def true_if_is_illumina(wildcards):
    return [{sample: is_illumina(None, sample)} for sample in get_samples_for_date(wildcards.date)]


wildcard_constraints:
    sample="[^/.]+",
    vartype="|".join(VARTYPES),
    clonality="subclonal|clonal",
    filter="|".join(
        list(map(re.escape, config["variant-calling"]["filters"])) + ["nofilter"]
    ),
    varrange="structural|small|homopolymer",
