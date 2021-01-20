from pathlib import Path


def get_samples():
    return pep.sample_table["sample_name"]


def get_fastqs(wildcards):
    return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]


def get_synonyms():
    return local(Path(workflow.snakefile).parent.parent / "resources/synonyms.txt")
