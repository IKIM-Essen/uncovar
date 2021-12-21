# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import re
import sys

sys.stderr = open(snakemake.log[0], "w")
# sys.stdout = open(snakemake.log[0], "a")

import gffutils
import numpy as np
import pandas as pd
import pysam


def phred_to_prob(phred):
    if phred is None:
        return 0
    return 10 ** (-phred / 10)


def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)


variants_df = pd.DataFrame()
lineage_df = pd.DataFrame()

# read generated variant file and extract all variants
with pysam.VariantFile(snakemake.input.variant_file, "rb") as infile:
    for record in infile:
        if "SIGNATURES" in record.info:
            signatures = record.info.get("SIGNATURES", ("#ERROR0",))[0]
            vaf = record.samples[0]["AF"][0]
            prob_not_present = phred_to_prob(
                record.info["PROB_ABSENT"][0]
            ) + phred_to_prob(record.info["PROB_ARTIFACT"][0])
            lineages = record.info["LINEAGES"]

            # generate df with all signatures + VAF and Prob_not_present from calculation
            variants_df = variants_df.append(
                {
                    "Signatures": signatures,
                    "VAF": vaf,
                    "Prob_not_present": prob_not_present,
                },
                ignore_index=True,
            )
            # generate df with lineage matrix for all signatures
            lineage_df = lineage_df.append(
                {"Signatures": signatures, **{lineage: "x" for lineage in lineages}},
                ignore_index=True,
            )

# count occurences of signatures (x) in lineage columns and get sorted list
lineage_dict = dict(lineage_df.count())
lineage_dict = dict(
    sorted(lineage_dict.items(), key=lambda item: item[1], reverse=True)
)
top5_lineages = list(lineage_dict.keys())

# only include variant names (index=0) + top 5 variants (index=1-6) and reorder
lineage_df.drop(labels=top5_lineages[7:], axis=1, inplace=True)
lineage_df = lineage_df[top5_lineages[:7]]

# aggregate both dataframes by summing up repeating rows for VAR (maximum=1) and multiply Prob_not_present
variants_df = (
    variants_df.groupby(["Signatures"])
    .agg(func={"VAF": lambda x: min(sum(x), 1.0), "Prob_not_present": np.prod}, axis=1)
    .reset_index()
)

# new column for 1-prob_not_present = prob_present
variants_df["Prob_present"] = 1.0 - variants_df["Prob_not_present"]
variants_df["Prob X VAF"] = variants_df["Prob_present"] * variants_df["VAF"]
lineage_df = lineage_df.groupby(["Signatures"]).agg("max").reset_index()

# calculate Jaccard coefficient for top 5 lineages and save row as df to append after sorting
jaccard_coefficient = {}
for lineage in range(1, len(top5_lineages[:7])):
    jaccard_coefficient[top5_lineages[lineage]] = (
        variants_df[
            variants_df["Signatures"].isin(
                lineage_df[lineage_df[top5_lineages[lineage]] == "x"]["Signatures"]
            )
        ]["Prob X VAF"].sum()
        / variants_df["Prob X VAF"].sum()
    )
jaccard_row = pd.DataFrame({"Signatures": "Jaccard", **jaccard_coefficient}, index=[0])

# merge variants dataframe and lineage dataframe
variants_df = variants_df.merge(lineage_df, left_on="Signatures", right_on="Signatures")

# add feature column for sorting
variants_df["Features"] = variants_df["Signatures"].str.extract(r"(.+)[:].+|\*")

# position of variant for sorting and change type
variants_df["Position"] = variants_df["Signatures"].str.extract(
    r"([0-9]+)([A-Z]+|\*)$"
)[0]
variants_df = variants_df.astype({"Position": "int64"})

# generate sorting list from .gff with correct order of features
gff = gffutils.create_db(snakemake.input.annotation, dbfn=":memory:")
gene_start = {gene["gene_name"][0]: gene.start for gene in gff.features_of_type("gene")}
sorter = [k[0] for k in sorted(gene_start.items(), key=lambda item: item[1])]
sorterIndex = dict(zip(sorter, range(len(sorter))))
variants_df["Features_Rank"] = variants_df["Features"].map(sorterIndex)

# define categories for sorting
variants_df.loc[
    (variants_df[top5_lineages[1]] == "x") & (variants_df["Prob_present"] >= 0.95),
    "Order",
] = 0
variants_df.loc[
    (variants_df[top5_lineages[1]] == "x") & (variants_df["Prob_present"] <= 0.05),
    "Order",
] = 1
variants_df.loc[
    (variants_df[top5_lineages[1]] == "x")
    & ((variants_df["Prob_present"] > 0.05) & (variants_df["Prob_present"] < 0.95)),
    "Order",
] = 2
variants_df.loc[
    (variants_df[top5_lineages[1]] != "x") & (variants_df["Prob_present"] <= 0.05),
    "Order",
] = 3
variants_df.loc[
    (variants_df[top5_lineages[1]] != "x") & (variants_df["Prob_present"] >= 0.95),
    "Order",
] = 4
variants_df.loc[
    (variants_df[top5_lineages[1]] != "x")
    & ((variants_df["Prob_present"] > 0.05) & (variants_df["Prob_present"] < 0.95)),
    "Order",
] = 5

# sort final DF
variants_df["Prob X VAF"].replace([0, 0.0], np.NaN, inplace=True)
variants_df.sort_values(
    by=["Order", "Features_Rank", "Position"],
    ascending=[True, True, True],
    na_position="last",
    inplace=True,
)

# concat row with Jaccard coefficient, drop unneccesary columns, sort with Jaccard coefficient, round
variants_df = pd.concat([jaccard_row, variants_df]).reset_index(drop=True)
variants_df = variants_df[["Signatures", "Prob_present", "VAF", *top5_lineages[1:7]]]
variants_df = variants_df.round({"Prob_present": 5, "VAF": 5})
variants_df.set_index("Signatures", inplace=True)
variants_df.sort_values(
    by="Jaccard", axis=1, na_position="first", ascending=False, inplace=True
)
# output variant_df
variants_df.to_csv(snakemake.output.variant_table, index=True, sep=",")

# TODO
# mutations of top lineage present, that are not present overall
# mutations of sample, not in top lineage
# mutations of top lineage, that are present over all
# all sorted by feature rank and postion per block
# for small probabilties fsum

# DONE
# sorted with top lineage, feature_rank, position
# Jaccard coefficient:
# multiply VAF with prob before
# prob of all top lineage mutations summed up divided by sum of all probs of all mutations
# keep VAF and prob
# function for min 1 for VAF
