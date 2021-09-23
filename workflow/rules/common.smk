from pathlib import Path
import pandas as pd
import re
import random


VARTYPES = ["SNV", "MNV", "INS", "DEL", "REP", "INV", "DUP"]

BENCHMARK_PREFIX = "benchmark-sample-"
NON_COV2_TEST_PREFIX = "non-cov2-"
MIXTURE_PREFIX = "mixture-sample-"
MIXTURE_PART_INDICATOR = "_MIX_"
MIXTURE_PERCENTAGE_INDICATOR = "_PERC_"
BENCHMARK_DATE_WILDCARD = "benchmarking"


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


def get_samples_before_date(wildcards):
    return list(
        pep.sample_table[pep.sample_table["run_id"] <= wildcards.date][
            "sample_name"
        ].values
    )


def get_dates_before_date(wildcards):
    return list(
        pep.sample_table[pep.sample_table["run_id"] <= wildcards.date]["run_id"].values
    )


def get_fastqs(wildcards):
    if wildcards.sample.startswith(BENCHMARK_PREFIX):
        # this is a simulated benchmark sample, do not look up FASTQs in the sample sheet
        accession = wildcards.sample[len(BENCHMARK_PREFIX) :]
        return expand(
            "resources/benchmarking/{accession}/reads.{read}.fastq.gz",
            accession=accession,
            read=[1, 2],
        )
    if wildcards.sample.startswith(NON_COV2_TEST_PREFIX):
        # this is for testing non-sars-cov2-genomes
        accession = wildcards.sample[len(NON_COV2_TEST_PREFIX) :]
        return expand(
            "resources/benchmarking/{accession}/reads.{read}.fastq.gz",
            accession=accession,
            read=[1, 2],
        )
    if wildcards.sample.startswith(MIXTURE_PREFIX):
        mixture = wildcards.sample[len(MIXTURE_PREFIX) :]
        return expand(
            "resources/mixtures/{mixtures}/reads.{read}.fastq.gz",
            mixtures=mixture,
            read=[1, 2],
        )
    # default case, look up FASTQs in the sample sheet
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
        with checkpoints.extract_strain_genomes_from_gisaid.get(
            date=wildcards.date
        ).output[0].open() as f:
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
    ):
        return expand(
            "results/{date}/trimmed/{sample}.{read}.fastq.gz",
            date=wildcards.date,
            read=[1, 2],
            sample=wildcards.sample,
        )

    # theses reads are used to generate the bam file for the BAMclipper
    elif wildcards.reference == config["adapters"]["amplicon-reference"]:
        return expand(
            "results/{date}/nonhuman-reads/{sample}.{read}.fastq.gz",
            date=wildcards.date,
            read=[1, 2],
            sample=wildcards.sample,
        )

    # aligments to other references (e.g. the covid reference genome),
    # are done with reads, which have undergone the quality control process
    else:
        return get_reads_after_qc(wildcards)


def get_reads_after_qc(wildcards, read="both"):

    if is_amplicon_data(wildcards.sample):
        pattern = expand(
            "results/{date}/clipped-reads/{sample}.{read}.fastq.gz",
            date=wildcards.date,
            read=[1, 2],
            sample=wildcards.sample,
        )
    else:
        pattern = expand(
            "results/{date}/nonhuman-reads/{sample}.{read}.fastq.gz",
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
    conf = config["RKI-quality-criteria"]
    if is_amplicon_data(wildcards.sample):
        return conf["min-depth-with-PCR-duplicates"]
    else:
        return conf["min-depth-without-PCR-duplicates"]


def return_assembler(sample):
    if is_amplicon_data(sample):
        return config["assemblers_used"]["amplicon"]
    else:
        return config["assemblers_used"]["shotgun"]


def get_contigs(wildcards):
    pattern = (
        "results/{date}/assembly/{sample}/{assembler}/{sample}.contigs.fasta".format(
            assembler=return_assembler(wildcards.sample), **wildcards
        ),
    )
    return pattern


def get_expanded_contigs(wildcards):
    sample = get_samples_for_date(wildcards.date)
    return [
        "results/{{date}}/assembly/{sample}/{assembler}/{sample}.contigs.fasta".format(
            sample=s, assembler=return_assembler(s)
        )
        for s in sample
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
    elif wildcards.stage == "masked":
        return "results/{date}/contigs/masked/{sample}.fasta"
    elif wildcards.stage == "pseudoassembly":
        return "results/{date}/contigs/pseudoassembled/{sample}.fasta"


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
    else:
        return max_reads


def get_strain(path_to_pangolin_call):
    pangolin_results = pd.read_csv(path_to_pangolin_call)
    return pangolin_results.loc[0]["lineage"]


def is_amplicon_data(sample):
    if (
        sample.startswith(BENCHMARK_PREFIX)
        or sample.startswith(NON_COV2_TEST_PREFIX)
        or sample.startswith(MIXTURE_PREFIX)
    ):
        # benchmark data, not amplicon based
        return False
    sample = pep.sample_table.loc[sample]
    try:
        return bool(int(sample["is_amplicon_data"]))
    except KeyError:
        return False


def get_samples_for_date_amplicon(date):
    return [s for s in get_samples_for_date(date) if is_amplicon_data(s)]


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
    if is_amplicon_data(wildcards.sample):
        # use clipped reads
        return "results/{date}/clipped-reads/{sample}.primerclipped.bam"
    # use trimmed reads
    amplicon_reference = config["adapters"]["amplicon-reference"]
    return "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam".format(
        ref=amplicon_reference
    )


def get_adapters(wildcards):
    if is_amplicon_data(wildcards.sample):
        return config["adapters"]["illumina-nimagen"]
    return config["adapters"]["illumina-revelo"]


def get_final_assemblies(wildcards):
    if wildcards.assembly_type == "masked-assembly":
        pattern = "results/{{date}}/contigs/masked/{sample}.fasta"
    elif wildcards.assembly_type == "pseudo-assembly":
        pattern = "results/{{date}}/contigs/pseudoassembled/{sample}.fasta"

    return expand(pattern, sample=get_samples_for_date(wildcards.date))


def get_final_assemblies_identity(wildcards):
    if wildcards.assembly_type == "masked-assembly":
        pattern = "results/{{date}}/quast/masked/{sample}/report.tsv"
    elif wildcards.assembly_type == "pseudo-assembly":
        pattern = "results/{{date}}/quast/pseudoassembly/{sample}/report.tsv"

    return expand(pattern, sample=get_samples_for_date(wildcards.date))


def get_assemblies_for_submission(agg_type):
    def inner(wildcards):
        if wildcards.date != BENCHMARK_DATE_WILDCARD:
            with checkpoints.rki_filter.get(
                date=wildcards.date, assembly_type="masked-assembly"
            ).output[0].open() as f:

                masked_samples = (
                    pd.read_csv(f, squeeze=True, header=None).astype(str).to_list()
                )

            with checkpoints.rki_filter.get(
                date=wildcards.date, assembly_type="pseudo-assembly"
            ).output[0].open() as f:
                pseudo_samples = (
                    pd.read_csv(f, squeeze=True, header=None).astype(str).to_list()
                )
        # for testing of pangolin don't create pseudo-assembly
        else:
            masked_samples = [wildcards.sample]

        pseudo_assembly_pattern = (
            "results/{{date}}/contigs/pseudoassembled/{sample}.fasta"
        )
        normal_assembly_pattern = "results/{{date}}/contigs/masked/{sample}.fasta"

        # get accepted samples for rki submission
        if agg_type == "accepted samples":
            accepted_assemblies = []

            for sample in set(masked_samples + pseudo_samples):
                if sample in masked_samples:
                    accepted_assemblies.append(
                        normal_assembly_pattern.format(sample=sample)
                    )
                else:
                    accepted_assemblies.append(
                        pseudo_assembly_pattern.format(sample=sample)
                    )
            return accepted_assemblies

        # for the pangolin call
        elif agg_type == "single sample":
            if wildcards.sample in masked_samples:
                return "results/{date}/contigs/polished/{sample}.fasta"
            elif wildcards.sample in pseudo_samples:
                return "results/{date}/contigs/pseudoassembled/{sample}.fasta"
            # for not accepted samples use the polished-contigs
            else:
                return "results/{date}/contigs/polished/{sample}.fasta"

        # for the qc report
        elif agg_type == "all samples":
            assembly_type_used = []
            for sample in get_samples_for_date(wildcards.date):
                if sample in masked_samples:
                    assembly_type_used.append(f"{sample},normal")
                elif sample in pseudo_samples:
                    assembly_type_used.append(f"{sample},pseudo")
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
            sample=get_samples_for_date(wildcards.date),
            **kwargs,
        )

    return inner


def expand_samples_for_date(paths, **kwargs):
    return expand_samples_by_func(paths, get_samples_for_date, **kwargs)


def expand_samples_for_date_amplicon(paths, **kwargs):
    return expand_samples_by_func(paths, get_samples_for_date_amplicon, **kwargs)


def get_vep_args(wildcards, input):
    return (
        "--vcf_info_field ANN --hgvsg --hgvs --synonyms {synonyms} "
        "--custom {input.problematic},,vcf,exact,0,"
    ).format(input=input, synonyms=get_resource("synonyms.txt"))


def expand_samples_by_func(paths, func, **kwargs):
    def inner(wildcards):
        return expand(
            paths,
            sample=get_samples_for_date(wildcards.date),
            **kwargs,
        )

    return inner


def get_samples_for_assembler_comparison(paths):
    return zip_expand(
        paths,
        get_dates(),
        get_samples(),
        config["assemblers_for_comparison"],
    )


def get_megahit_preset(wildcards):
    if wildcards.preset == "std":
        return ""
    else:
        return f"--preset {wildcards.preset}"


wildcard_constraints:
    sample="[^/.]+",
    vartype="|".join(VARTYPES),
    clonality="subclonal|clonal",
    filter="|".join(
        list(map(re.escape, config["variant-calling"]["filters"])) + ["nofilter"]
    ),
    varrange="structural|small",
