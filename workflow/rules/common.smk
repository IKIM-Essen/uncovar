def get_samples():
    return pep.sample_table["sample_name"]


def get_fastqs(wildcards):
    return pep.sample_table.loc[wildcards.sample][["fq1", "fq2"]]
