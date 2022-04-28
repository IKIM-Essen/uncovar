# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import re
import sys

sys.stderr = open(snakemake.log[0], "w")

import gffutils
import numpy as np
import pandas as pd
import pysam


def phred_to_prob(phred):
    if phred is None:
        return pd.NA
    return 10 ** (-phred / 10)


# np.prod returns 1's as values for a pd series with NaN's. A list would return NaN's
def prod_prob_not_present(probs):
    if pd.isna(probs).any():
        return pd.NA
    else:
        return np.prod(probs)


def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)


def add_number_suffix(number):
    number = str(number)

    if number.endswith("1") and number != "11":
        return f"{number}st"
    elif number.endswith("2") and number != "12":
        return f"{number}nd"
    elif number.endswith("3") and number != "13":
        return f"{number}rd"
    else:
        return f"{number}th"


def rename_enumeration(list_length):
    return [add_number_suffix(x) for x in range(1, list_length + 1)]


variants_df = pd.DataFrame()
lineage_df = pd.DataFrame()

# read generated variant file and extract all variants
with pysam.VariantFile(snakemake.input.variant_file, "rb") as infile:
    for record in infile:
        if "SIGNATURES" in record.info:
            signatures = record.info.get("SIGNATURES", ("#ERROR0",))
            vaf = record.samples[0]["AF"][0]
            dp = record.samples[0]["DP"]
            prob_not_present = phred_to_prob(
                record.info["PROB_ABSENT"][0]
            ) + phred_to_prob(record.info["PROB_ARTIFACT"][0])
            if pd.isna(prob_not_present):
                vaf = pd.NA
            lineages = record.info["LINEAGES"]
            for signature in signatures:
                # generate df with all signatures + VAF and Prob_not_present from calculation
                # variants_df = variants_df.append(
                #    {
                #        "Mutations": signature,
                #        "Frequency": vaf,
                #        "ReadDepth": dp,
                #        "Prob_not_present": prob_not_present,
                #    },
                #    ignore_index=True,
                # )
                variants_df_append = {
                    "Mutations": signature,
                    "Frequency": vaf,
                    "ReadDepth": dp,
                    "Prob_not_present": prob_not_present,
                }
                variants_df = pd.concat(
                    [variants_df, variants_df_append], ignore_index=True
                )
                # generate df with lineage matrix for all signatures
                # lineage_df = lineage_df.append(
                #    {
                #        "Mutations": signature,
                #        **{lineage.replace(".", " "): "x" for lineage in lineages},
                #    },
                #    ignore_index=True,
                # )
                lineage_df_append = {
                    "Mutations": signature,
                    **{lineage.replace(".", " "): "x" for lineage in lineages},
                }
                lineage_df = pd.concat(
                    [lineage_df, lineage_df_append], ignore_index=True
                )

# aggregate both dataframes by summing up repeating rows for VAR (maximum=1) and multiply Prob_not_present
variants_df = variants_df.groupby(["Mutations"]).agg(
    func={
        "Frequency": lambda x: min(sum(x), 1.0),
        "Prob_not_present": prod_prob_not_present,
        "ReadDepth": np.min,
    },
    axis=1,
)

# new column for 1-prob_not_present = prob_present
variants_df["Probability"] = 1.0 - variants_df["Prob_not_present"]
variants_df["Prob X VAF"] = variants_df["Probability"] * variants_df["Frequency"]
lineage_df = lineage_df.drop_duplicates()

# merge duplicated mutations
lineage_df = lineage_df.fillna(0)
lineage_df = lineage_df.replace({"x": 1})
lineage_df = (
    lineage_df.groupby(["Mutations"])
    .agg(func={column: np.max for column in lineage_df.columns})
    .reset_index(drop=True)
)
lineage_df = lineage_df.replace({1: "x", 0: ""})

# calculate Jaccard coefficient for each lineage
# iterate over lineages in columns (mutations as index)
lineage_df.set_index("Mutations", inplace=True)
jaccard_coefficient = {}
for lineage in lineage_df.columns:
    lineage_defining_variants = variants_df.index.isin(
        lineage_df.index[lineage_df[lineage] == "x"]
    )
    lineage_defining_non_variants = ~lineage_defining_variants
    jaccard_coefficient[lineage] = round(
        (
            variants_df[lineage_defining_variants]["Prob X VAF"].sum()
            + variants_df[lineage_defining_non_variants]["Prob_not_present"].sum()
        )
        / len(variants_df),
        3,
    )

jaccard_row = pd.DataFrame(
    {"Mutations": "Similarity", **jaccard_coefficient}, index=[0]
)
# merge variants dataframe and lineage dataframe
variants_df = variants_df.merge(lineage_df, left_index=True, right_index=True)

# add feature column for sorting
variants_df["Features"] = variants_df.index.to_series().str.extract(r"(.+)[:].+|\*")

# position of variant for sorting and change type
variants_df["Position"] = variants_df.index.to_series().str.extract(
    r"(.*:?[A-Z]+|\*$|-)([0-9]+)([A-Z]+$|\*$|-)$"
)[1]
variants_df = variants_df.astype({"Position": "int64"})

# generate sorting list from .gff with correct order of features
gff = gffutils.create_db(snakemake.input.annotation, dbfn=":memory:")
gene_start = {gene["gene_name"][0]: gene.start for gene in gff.features_of_type("gene")}
sorter = [k[0] for k in sorted(gene_start.items(), key=lambda item: item[1])]
sorterIndex = dict(zip(sorter, range(len(sorter))))
variants_df["Features_Rank"] = variants_df["Features"].map(sorterIndex)

# row for lineage name after renaming columns (column names can't be formatted)
lineages_row_df = pd.DataFrame(
    {
        "Mutations": "Lineage",
        **{x: x for x in list(lineage_df.columns) if x != "Mutations"},
    },
    index=[0],
)

# concat row with Jaccard coefficient, drop unneccesary columns, sort with Jaccard coefficient, round
variants_df.reset_index(inplace=True)
variants_df = pd.concat([jaccard_row, variants_df])
variants_df = pd.concat([lineages_row_df, variants_df])


variants_df = variants_df.round({"Probability": 2, "Frequency": 2})
variants_df.set_index("Mutations", inplace=True)
variants_df.sort_values(
    by="Similarity", axis=1, na_position="first", ascending=False, inplace=True
)
# rename hits ascending
variants_df.rename(
    columns={
        x: y
        for x, y in zip(
            list(variants_df.columns)[8:],
            rename_enumeration(len(list(variants_df.columns)[8:])),
        )
    },
    errors="raise",
    inplace=True,
)
# sort final DF
variants_df.loc[
    variants_df["1st"] == "x",
    "Order",
] = 1
variants_df.loc[
    variants_df["1st"] != "x",
    "Order",
] = 2
variants_df.at["Similarity", "Order"] = 0
variants_df.at["Lineage", "Order"] = 0
variants_df["Prob X VAF"].replace([0, 0.0], np.NaN, inplace=True)
variants_df.sort_values(
    by=["Order", "Features_Rank", "Position"],
    ascending=[True, True, True],
    na_position="last",
    inplace=True,
)
# drop unwanted columns
variants_df.drop(
    columns=[
        "Prob_not_present",
        "Prob X VAF",
        "Features",
        "Position",
        "Features_Rank",
        "Order",
    ],
    inplace=True,
)
all_columns = variants_df.columns
first_columns = ["Probability", "Frequency", "ReadDepth"]
rest_columns = [item for item in all_columns if item not in first_columns]
variants_df = variants_df[[*first_columns, *rest_columns]]

# drop other lineages, top 10 only
variants_df.drop(variants_df.columns[13:], axis=1, inplace=True)

# output variant_df
variants_df.to_csv(snakemake.output.variant_table, index=True, sep=",")
