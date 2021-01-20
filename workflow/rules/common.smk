from pathlib import Path


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


def get_report_args(wildcards, files):
    return expand(
        "{sample}={file}", zip, sample=get_report_samples(wildcards), file=files
    )


def get_report_samples(wildcards):
    return get_samples() if wildcards.target == "all" else [wildcards.target]


def get_merge_calls_input(suffix):
    def inner(wildcards):
        return expand(
            "results/filtered-calls/{{sample}}.{vartype}.bcf.{suffix}",
            suffix=suffix,
            vartype=["SNV", "MNV", "INS", "DEL", "REP"],
        )
    return inner



wildcard_constraints:
    sample="|".join(get_samples())