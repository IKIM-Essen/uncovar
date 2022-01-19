# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import numpy as np


def get_test_cases_variant_calls(technology, suffix="", get="path"):
    """Returns bcf file paths used for generating varlociraptor test cases."""

    def inner(wildcards):
        sample_table = pep.sample_table.copy()
        sample_table = sample_table.loc[
            (sample_table["technology"] == technology)
            & (sample_table["test_case"] == wildcards.test_case)
        ]
        sample_table.sort_values(by=["test_case", "technology"], inplace=True)

        assert (
            len(sample_table) == 1
        ), f'Too many sampels are defined with technology "{technology}" for test case {wildcards.test_case}.'

        bcf_path_high = "results/{date}/filtered-calls/ref~main/{sample}.subclonal.high+moderate-impact.bcf{suffix}"
        bcf_path_low = "results/{date}/filtered-calls/ref~main/{sample}.subclonal.low-impact.bcf{suffix}"

        high_impact = expand(
            bcf_path_high,
            zip,
            date=sample_table["date"],
            sample=sample_table["sample_name"],
            suffix=suffix,
        )[0]
        low_impact = expand(
            bcf_path_low,
            zip,
            date=sample_table["date"],
            sample=sample_table["sample_name"],
            suffix=suffix,
        )[0]

        if get == "path":
            return [high_impact, low_impact]
        if get == "date":
            return sample_table["date"].to_list()[0]
        if get == "sample":
            return sample_table["sample_name"].to_list()[0]

    return inner


def get_test_cases_data(technology, typ):
    sample_table = pep.sample_table.dropna(subset=["test_case"])
    filter_by = sample_table["technology"] == technology
    if typ == "sample":
        return sample_table.loc[filter_by].index.to_list()
    if typ == "date":
        return sample_table.loc[filter_by]["date"].to_list()

    raise TypeError(f"Technology {technology} or type {typ} recognzied")


def get_all_test_cases_names(wildcards):
    try:
        return pep.sample_table["test_case"].dropna().unique()
    except KeyError:
        raise TypeError("Column test_case not found in sample sheet.")


def get_aggregated_test_case_variants(return_typ):
    def inner(wildcards):
        variants = []
        test_case_paths = []
        vcf_paths = []
        csi_paths = []
        poses = []

        with checkpoints.get_test_case_variant_paths.get().output.paths.open() as f:
            for line in f.read().splitlines():
                pos, variant, test_case_path, vcf_path = line.split("\t")
                poses.append(pos)
                variants.append(variant)
                test_case_paths.append(test_case_path)
                vcf_paths.append(vcf_path)
                csi_paths.append(f"{vcf_path}.csi")

        if return_typ == "bcf":
            return vcf_paths
        if return_typ == "csi":
            return csi_paths
        elif return_typ == "variants":
            return variants
        elif return_typ == "poses":
            return poses
        elif return_typ == "test-case-paths":
            return test_case_paths

        raise TypeError(f"{return_typ} not recognized.")

    return inner


def get_test_cases(wildcards):
    with checkpoints.check_presence_of_test_case_variant_in_call.get().output[
        0
    ].open() as f:
        return f.read().splitlines()


def get_barcode_path(wildcards):
    return pep.sample_table.loc[wildcards.sample]["barcode"]
