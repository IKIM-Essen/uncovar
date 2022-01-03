# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from snakemake.io import expand

test_case_path = "{pos}\t{variant}\tresults/{date}/call-test-cases/ref~main/{sample}.{{varrange}}.chrom~{chrom}.pos~{pos}.bcf\tresults/{date}/candidate-calls/ref~main/{sample}.{{varrange}}.bcf"

variants = []
paths = []

for test_case in snakemake.input:
    variants_to_test = pd.read_csv(test_case, sep="\t")
    variants.append(variants_to_test)

variants = pd.concat(variants)

if len(variants) > 0:
    variants.loc[
        variants["vaf_illumina"] < variants["vaf_ont"], "to_test"
    ] = snakemake.params.illumina
    variants.loc[
        variants["vaf_illumina"] > variants["vaf_ont"], "to_test"
    ] = snakemake.params.ont

    sample_table = pd.DataFrame(snakemake.params.sample_table).dropna(
        subset=["test_case"]
    )
    sample_table = sample_table[["test_case", "technology", "sample_name", "date"]]

    variants["test_case"] = variants["test_case"].astype(str)
    variants = pd.merge(
        variants,
        sample_table,
        how="left",
        left_on=["test_case", "to_test"],
        right_on=["test_case", "technology"],
    )

    illumina_variants = variants.loc[variants["to_test"] == snakemake.params.illumina]
    illumina_paths = expand(
        test_case_path,
        zip,
        date=illumina_variants["date"],
        sample=illumina_variants["sample_name"],
        chrom=illumina_variants["chrom"],
        pos=illumina_variants["pos"],
        variant=illumina_variants["variant"],
    )
    illumina_paths = expand(illumina_paths, varrange=snakemake.params.illumina_varrange)
    paths.append(illumina_paths)

    ont_variants = variants.loc[variants["to_test"] == snakemake.params.ont]
    ont_variants = expand(
        test_case_path,
        zip,
        date=ont_variants["date"],
        sample=ont_variants["sample_name"],
        chrom=ont_variants["chrom"],
        pos=ont_variants["pos"],
        variant=ont_variants["variant"],
    )
    ont_variants = expand(ont_variants, varrange=snakemake.params.ont_varrange)
    paths.append(ont_variants)

    paths = sum(paths, [])

with open(snakemake.output.paths, "w") as f:
    for path in paths:
        print(f"{path}", file=f)

if len(variants) > 0:
    variants.drop(columns=["technology"]).to_csv(
        snakemake.output.overview, sep="\t", index=False
    )
else:
    pd.DataFrame().to_csv(snakemake.output.overview, sep="\t", index=False)
