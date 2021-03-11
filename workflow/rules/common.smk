from pathlib import Path
import pandas as pd
import re
import random


VARTYPES = ["SNV", "MNV", "INS", "DEL", "REP"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_dates():
    return list(pep.sample_table["run_id"].values)


def get_samples_for_date(date, filtered=False):
    # select samples with given date
    df = pep.sample_table
    df = df[df["run_id"] == date]

    samples_of_run = list(df["sample_name"].values)

    # filter
    if filtered:
        with checkpoints.rki_filter.get(date=date).output[0].open() as f:

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
    sorted_list = list(pep.sample_table["run_id"].unique())
    sorted_list.sort()
    return sorted_list


def get_latest_run_date():
    return pep.sample_table["run_id"].max()


def get_fastqs(wildcards, benchmark_prefix="benchmark-sample-"):
    if wildcards.sample.startswith(benchmark_prefix):
        # this is a simulated benchmark sample, do not look up FASTQs in the sample sheet
        accession = wildcards.sample[len(benchmark_prefix) :]
        return expand(
            "resources/benchmarking/{accession}/reads.{read}.fastq.gz",
            accession=accession,
            read=[1, 2],
        )
    if wildcards.sample.startswith("non-cov2-"):
        # this is for testing non-sars-cov2-genomes
        accession = wildcards.sample[len("non-cov2-") :]
        return expand(
            "resources/test-cases/{accession}/reads.{read}.fastq.gz",
            accession=accession,
            read=[1, 2],
        )
    if wildcards.sample.startswith("mixture-sample-"):
        print(wildcards)
        accessions = wildcards.sample[len("mixture-sample-") :]
        print(accessions)
        return expand(
            "resources/mixtures/{accession}/reads.{read}.fastq.gz",
            accession=accession,
            read=[1, 2],
        )
    # default case, look up FASTQs in the sample sheet
    return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]


def get_resource(name):
    return str((Path(workflow.snakefile).parent.parent / "resources") / name)


def get_report_input(pattern):
    def inner(wildcards):
        return expand(pattern, sample=get_report_samples(wildcards), **wildcards)

    return inner


def get_report_bcfs(wildcards, input):
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
        accessions = config.get("benchmark-genomes", [])
        if not accessions:
            accessions = pd.read_csv(f, squeeze=True)
        try:
            accessions = accessions[: config["limit-strain-genomes"]]
        except KeyError:
            # take all strain genomes
            pass
        return accessions


def get_non_cov2_accessions():
    accessions = config.get("non_cov2_genomes", [])
    return accessions


def get_strain_genomes(wildcards):
    # Case 1: take custom genomes from gisaid
    custom_genomes = config["strain-calling"]["use-gisaid"]
    if custom_genomes:
        with checkpoints.extract_strain_genomes_from_gisaid.get().output[0].open() as f:
            strain_genomes = pd.read_csv(f, squeeze=True).to_list()
            strain_genomes.append("resources/genomes/main.fasta")
            return expand("{strains}", strains=strain_genomes)

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
            "results/benchmarking/assembly/{accession}.bam"
            if bams
            else "resources/genomes/{accession}.fasta"
        )
        return expand(pattern, accession=accessions,)

    return inner


def get_non_cov2_calls(from_caller="pangolin"):
    accessions = get_non_cov2_accessions()
    pattern = (
        "results/test-cases/tables/strain-calls/non-cov2-{accession}.strains.pangolin.csv"
        if from_caller == "pangolin"
        else "results/test-cases/tables/strain-calls/non-cov2-{accession}.strains.kallisto.tsv"
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
            return "results/{date}/polished-contigs/{sample}.fasta".format(
                sample=wildcards.reference.replace("polished-", ""), **wildcards
            )
        else:
            # return assembly result
            return "results/{date}/ordered-contigs/{reference}.fasta{suffix}".format(
                suffix=suffix, **wildcards
            )

    return inner


def get_reads(wildcards):
    if (
        wildcards.reference == "human"
        or wildcards.reference == "main+human"
        or wildcards.reference.startswith("polished-")
    ):
        # alignment against the human reference genome must be done with trimmed reads,
        # since this alignment is used to generate the ordered, non human contigs
        return expand(
            "results/{date}/trimmed/{sample}.{read}.fastq.gz",
            date=wildcards.date,
            read=[1, 2],
            sample=wildcards.sample,
        )
    else:
        # other reference (e.g. the covid reference genome, are done with contigs) that
        # do not contain human contaminations
        return expand(
            "results/{date}/nonhuman-reads/{sample}.{read}.fastq.gz",
            date=wildcards.date,
            read=[1, 2],
            sample=wildcards.sample,
        )


def get_bwa_index(wildcards):
    if wildcards.reference == "human" or wildcards.reference == "main+human":
        return rules.bwa_large_index.output
    else:
        return rules.bwa_index.output


def get_target_events(wildcards):
    if wildcards.reference == "main" or wildcards.clonality != "clonal":
        # calling variants against the wuhan reference or we are explicitly interested in subclonal as well
        return "SUBCLONAL CLONAL"
    else:
        # only keep clonal variants
        return "CLONAL"


def get_filter_odds_input(wildcards):
    if wildcards.reference == "main":
        return "results/{date}/filtered-calls/ref~{reference}/{sample}.{filter}.bcf"
    else:
        return "results/{date}/calls/ref~{reference}/{sample}.bcf"


def get_vembrane_expression(wildcards):
    return config["variant-calling"]["filters"][wildcards.filter]


def zip_expand(expand_string, zip_wildcard_1, zip_wildcard_2, expand_wildcard):
    """
    Zip by two wildcards and the expand the zip over another wildcard.
    expand_string must contain {zip1}, {zip2} and {exp}.
    """

    return sum(
        [
            expand(ele, exp=expand_wildcard)
            for ele in expand(
                expand_string, zip, zip1=zip_wildcard_1, zip2=zip_wildcard_2,
            )
        ],
        [],
    )


def get_quast_fastas(wildcards):
    if wildcards.stage == "unpolished":
        return "results/{date}/assembly/{sample}/{sample}.contigs.fa"
    elif wildcards.stage == "polished":
        return "results/{date}/polished-contigs/{sample}.fasta"
    elif wildcards.stage == "masked":
        return "results/{date}/contigs-masked/{sample}.fasta"


def get_random_strain():
    with checkpoints.extract_strain_genomes_from_gisaid.get().output[0].open() as f:
        lines = f.read().splitlines()
        rnd_strain_path = random.choice(lines)
        strain = rnd_strain_path.split('/')[-1].replace(".fasta","")
        return strain


def get_mixture_results(wildcards):
    no_mixtures = config["mixtures"]["no_mixtures"]
    no_strains = config["mixtures"]["no_strains"]
    mixture_list = []

    for mix in range(no_mixtures):

        fractions = [random.randint(1, 100) for _ in range(no_strains)]
        s = sum(fractions)
        fractions = [ round(i/s *100) for i in fractions ]

        s = sum(fractions)
        if s!=100:
            fractions[-1] += (100 - s)

        mixture = ""
        for frac in fractions:
            strain = get_random_strain()
            mixture += f"#{strain}={frac}"
        
        mixture_list.append(mixture.replace(".", "-"))

    return  expand("results/benchmarking/tables/strain-calls/mixture-sample-{mixtures}.strains.kallisto.tsv", mixtures = mixture_list)


wildcard_constraints:
    sample="[^/.]+",
    vartypes="|".join(VARTYPES),
    clonality="subclonal|clonal",
    filter="|".join(
        list(map(re.escape, config["variant-calling"]["filters"])) + ["nofilter"]
    ),
