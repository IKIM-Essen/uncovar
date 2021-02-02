from pathlib import Path
import pandas as pd


VARTYPES = ["SNV", "MNV", "INS", "DEL", "REP"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_fastqs(wildcards, benchmark_prefix="benchmark-sample-"):
    if wildcards.sample.startswith(benchmark_prefix):
        # this is a simulated benchmark sample, do not look up FASTQs in the sample sheet
        accession = wildcards.sample[len(benchmark_prefix) :]
        return expand(
            "resources/benchmarking/{accession}/reads.{read}.fastq.gz",
            accession=accession,
            read=[1, 2],
        )
    # default case, look up FASTQs in the sample sheet
    return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]


def get_resource(name):
    return str((Path(workflow.snakefile).parent.parent / "resources") / name)


def get_report_input(pattern):
    def inner(wildcards):
        return expand(pattern, sample=get_report_samples(wildcards))

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
    return get_samples() if wildcards.target == "all" else [wildcards.target]


def get_merge_calls_input(suffix):
    def inner(wildcards):
        return expand(
            "results/filtered-calls/ref~{{reference}}/{{sample}}.{{clonality}}.{vartype}.fdr-controlled{suffix}",
            suffix=suffix,
            vartype=VARTYPES,
        )

    return inner


def get_strain_accessions(wildcards):
    with checkpoints.get_strain_accessions.get().output[0].open() as f:
        accessions = pd.read_csv(f, squeeze=True)
        try:
            accessions = accessions[: config["limit-strain-genomes"]]
        except KeyError:
            # take all strain genomes
            pass
        return accessions


def get_strain_genomes(wildcards):
    # Case 1: take custom genomes from config
    custom_genomes = config["strain-calling"].get("genomes", [])
    if custom_genomes:
        return custom_genomes

    # Case 2: take genomes from genbank
    accessions = get_strain_accessions(wildcards)
    return expand("resources/genomes/{accession}.fasta", accession=accessions)


def get_strain_signatures(wildcards):
    return expand(
        "resources/genomes/{accession}.sig", accession=get_strain_accessions(wildcards)
    )


def get_benchmark_results(wildcards):
    accessions = get_strain_accessions(wildcards)
    return expand(
        "results/tables/strain-calls/benchmark-sample-{accession}.strains.kallisto.tsv",
        accession=accessions,
    )


def get_aligned_contigs(wildcards):
    return expand("results/ordered_contigs/{sample}/{sample}.bam", sample=get_samples())


def get_assembly_contigs(wildcards):
    return expand("results/assembly/{sample}/final.contigs.fa", sample=get_samples())


def get_reference(suffix=""):
    def inner(wildcards):
        if wildcards.reference == "main":
            # return reference genome
            return "resources/genomes/main.fasta" + suffix
        else:
            # return assembly result
            return "results/ordered_contigs/{reference}.fasta.{suffix}".format(suffix=suffix, **wildcards)
    return inner


def get_target_events(wildcards):
    if wildcards.reference == "main" or wildcards.clonality != "clonal":
        # calling variants against the wuhan reference or we are explicitly interested in subclonal as well
        return "SUBCLONAL CLONAL"
    else:
        # only keep clonal variants
        return "CLONAL"


def get_filter_odds_input(wildcards):
    if wildcards.reference == "main":
        return "results/annotated-calls/ref~main/{sample}.bcf"
    else:
        return "results/calls/ref~main/{sample}.bcf"


wildcard_constraints:
    sample="[^/.]+",
    vartypes="|".join(VARTYPES),
    clonality="subclonal|clonal"
