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

gff = gffutils.create_db(snakemake.input.annotation, dbfn=":memory:")
gene_start = {gene["gene_name"][0]: gene.start for gene in gff.features_of_type("gene")}
sorter = [k[0] for k in sorted(gene_start.items(), key=lambda item: item[1])]

with pysam.VariantFile(snakemake.input[0], "rb") as infile:
    for record in infile:
        if "SIGNATURES" in record.info:
            signatures = record.info.get("SIGNATURES", ("#ERROR0",))[0]
            vaf = record.samples[0]["AF"][0]
            prob_not_present = phred_to_prob(
                record.info["PROB_ABSENT"][0]
            ) + phred_to_prob(record.info["PROB_ARTIFACT"][0])
            lineages = record.info["LINEAGES"]
            lineage_dict = {}
            # print(lineages)
            for item in lineages:
                lineage_dict[item] = "x"
            # print(lineage_dict)
            variants_df = variants_df.append(
                {
                    "Signatures": signatures,
                    "VAF": vaf,
                    "Prob_not_present": prob_not_present,
                },
                ignore_index=True,
            )
            lineage_df = lineage_df.append(
                {"Signatures": signatures, **lineage_dict}, ignore_index=True
            )

# count occurences of x in lineage columns and get sorted list
lineage_dict = dict(lineage_df.count())
lineage_dict = dict(
    sorted(lineage_dict.items(), key=lambda item: item[1], reverse=True)
)
top5_lineages = list(lineage_dict.keys())

# only include variant names + top 5 variants and reorder
lineage_df.drop(labels=top5_lineages[7:], axis=1, inplace=True)
lineage_df = lineage_df[top5_lineages[:7]]
lineage_df.to_csv(snakemake.output.lineage_df)


# aggregate both dataframes by summing up repeating rows for variants
variants_df = (
    variants_df.groupby(["Signatures"])
    .agg(func={"VAF": lambda x: min(sum(x), 1.0), "Prob_not_present": np.prod}, axis=1)
    .reset_index()
)
# new column for 1-prob_not_present = prob_present
variants_df["Prob_present"] = 1.0 - variants_df["Prob_not_present"]
variants_df["Prob X VAF"] = variants_df["Prob_present"] * variants_df["VAF"]
lineage_df = lineage_df.groupby(["Signatures"]).agg("max").reset_index()

# calculate Jaccard coefficient for top 5 lineages
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

# print(variants_df[variants_df["Signatures"] == top_df["Signatures"].values()])


# merge variants dataframe and lineage dataframe
variants_df = variants_df.merge(lineage_df, left_on="Signatures", right_on="Signatures")

# add feature column for sorting
variants_df["Features"] = variants_df["Signatures"].str.extract(r"(.+)[:].+|\*")

print(variants_df)
# position of variant for sorting and change type
variants_df["Position"] = variants_df["Signatures"].str.extract(
    r"([0-9]+)([A-Z]+|\*)$"
)[0]
variants_df = variants_df.astype({"Position": "int64"})

# generate sorting list with correct range of features
sorterIndex = dict(zip(sorter, range(len(sorter))))
variants_df["Features_Rank"] = variants_df["Features"].map(sorterIndex)

# replace zeros with empty string and sort final DF
# print(variants_df)
variants_df["Prob X VAF"].replace([0, 0.0], np.NaN, inplace=True)
# print(variants_df)
variants_df.sort_values(
    by=["Prob X VAF", top5_lineages[1], "Features_Rank", "Position"],
    ascending=[False, True, True, True],
    na_position="last",
    inplace=True,
)
variants_df = pd.concat([jaccard_row, variants_df]).reset_index(drop=True)
variants_df = variants_df[["Signatures", "Prob_present", "VAF", *top5_lineages[1:7]]]
variants_df.set_index("Signatures", inplace=True)
variants_df.sort_values(
    by="Jaccard", axis=1, na_position="first", ascending=False, inplace=True
)
variants_df.to_csv(snakemake.output[0], index=True, sep=",")

# mutations of top lineage present, that are not present overall
# mutations of sample, not in top lineage
# mutations of top lineage, that are present over all
# all sorted by feature rank and postion per block

# Jaccard coefficient:
# multiply VAF with prob before
# prob of all top lineage mutations summed up divided by sum of all probs of all mutations
# for small probabilties fsum
# keep VAF and prob
# function for min 1 for VAF
# function for
