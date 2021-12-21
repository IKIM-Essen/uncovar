import numpy as np

def get_test_cases_variant_calls(technology, suffix="", get="path"):
    """Returns bcf file paths used for generating varlociraptor test cases."""
    def inner(wildcards):
        sample_table = pep.sample_table.copy()
        sample_table =sample_table.loc[(sample_table["technology"] == technology) & (sample_table["test_case"] == wildcards.test_case)]
        sample_table.sort_values(by=["test_case", "technology"], inplace=True)

        assert len(sample_table) == 1, f"Too many sampels are defined with technology \"{technology}\" for test case {wildcards.test_case}."

        bcf_path="results/{date}/filtered-calls/ref~main/{sample}.subclonal.high+moderate-impact.bcf{suffix}"

        if get=="path":
            return expand(bcf_path, zip, date=sample_table["date"], sample=sample_table["sample_name"], suffix = suffix)[0]
        if get=="date":
            return sample_table["date"].to_list()[0]
        if get=="sample":
            return sample_table["sample_name"].to_list()[0]

    return inner


def get_test_cases_data(technology, typ):
    sample_table = pep.sample_table
    filter_by = (sample_table["technology"] == technology)
    if typ == "sample":
        return sample_table.loc[filter_by].index.to_list()[0]
    if typ == "date":
        return sample_table.loc[filter_by]["date"].to_list()[0]

    raise TypeError(f"Technology {technology} or type {typ} recognzied")


def get_all_test_cases_names(wildcards):
    try:
        return pep.sample_table["test_case"].dropna().unique()
    except KeyError:
        raise TypeError("Column test_case not found in sample sheet.")


def get_test_cases(wildcards):
    with checkpoints.aggregate_test_case_variants.get().output.paths.open() as f:
        return f.read().splitlines()
