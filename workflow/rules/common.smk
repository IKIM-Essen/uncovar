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
            "results/filtered-calls/{{sample}}.{vartype}{suffix}",
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
    accessions = get_strain_accessions(wildcards)
    return expand("resources/genomes/{accession}.fasta", accession=accessions)


def get_strain_signatures(wildcards):
    return expand(
        "resources/genomes/{accession}.sig", accession=get_strain_accessions(wildcards)
    )


def get_whitelist_strain_accessions(wildcards):
    with checkpoints.cat_genomes.get().output[0].open() as f:
        whitelist = config["whitelisting"]["whitelist-lineage"]

        # downloadable from https://www.epicov.org/epi3 / EpiCoV™ / Downloads / nextmeta. Must first register
        # TODO automate download of metadata file
        metadata = pd.read_csv(
            "resources/metadata.tsv", delimiter="\t", low_memory=False
        )
        metadata.replace("?", float("NaN"), inplace=True)
        metadata = metadata[metadata.pangolin_lineage.isin(whitelist)]
        metadata["cut_strain"] = metadata.strain.str.extract(r"(\/(.*)\/\d{4}$)")[1]
        metadata.set_index("strain", drop=True, inplace=True)

        whitelisted_accessions = []
        for genome in expand(
            "resources/genomes/{accession}.fasta", accession=wildcards
        ):
            with open(genome, "r") as f:
                first_line = f.readline()
                if any(
                    strain in first_line for strain in metadata.cut_strain.to_list()
                ):
                    whitelisted_accessions.append(first_line)
        whitelisted_accessions = (
            pd.Series(whitelisted_accessions)
            .str.extract(r">([^\s]+)\.")
            .dropna()[0]
            .unique()
        )
        return list(whitelisted_accessions)


def get_whitelisted_strain_genomes(wildcards):
    accessions = get_strain_accessions(wildcards)
    accessions = get_whitelist_strain_accessions(accessions)
    return expand("resources/genomes/{accession}.fasta", accession=accessions)


def get_concatenated_genomes(wildcards):
    if config["whitelisting"]["use-whitelist"]:
        return "resources/strain-genomes-whitelisted.fasta"
    else:
        return "resources/strain-genomes.fasta"


def get_benchmark_results(wildcards):
    accessions = get_strain_accessions(wildcards)
    return expand(
        "results/tables/strain-calls/benchmark-sample-{accession}.strains.tsv",
        accession=accessions,
    )


wildcard_constraints:
    sample="[^/.]+",
    vartypes="|".join(VARTYPES),
