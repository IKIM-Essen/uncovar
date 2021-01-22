from pathlib import Path

VARTYPES = ["SNV", "MNV", "INS", "DEL", "REP"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_fastqs(wildcards):
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

def get_strain_accessions_from_txt(file):
    with open(file, 'r') as f:
        strain_accessions = f.read().splitlines()
    strain_accessions.remove('id')
    return list(strain_accessions)[:100]

wildcard_constraints:
    sample="|".join(get_samples()),
    vartypes="|".join(VARTYPES),
