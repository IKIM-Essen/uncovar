# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import re
import sys

sys.stderr = open(snakemake.log[0], "w")
# sys.stdout = open(snakemake.log[0], "a")

import gffutils
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
            prob_clonal = phred_to_prob(record.info["PROB_CLONAL"][0])
            prob_subclonal_min = phred_to_prob(record.info["PROB_SUBCLONAL_MINOR"][0])
            prob_subclonal_maj = phred_to_prob(record.info["PROB_SUBCLONAL_MAJOR"][0])
            prob_subclonal_high = phred_to_prob(record.info["PROB_SUBCLONAL_HIGH"][0])
            prob_low = phred_to_prob(record.info["PROB_LOW"][0])
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
                    "Prob_clonal": prob_clonal,
                    # [x for x in lineages]: ["x" f],
                    # "Prob_subclonal_min": prob_subclonal_min,
                    # "Prob_subclonal_maj": prob_subclonal_maj,
                    # "Prob_subclonal_high": prob_subclonal_high,
                    # "Prob_low": prob_low,
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
keys = list(lineage_dict.keys())

# only include variant names + top 5 variants and reorder
lineage_df.drop(labels=keys[7:], axis=1, inplace=True)
lineage_df = lineage_df[keys[:7]]
lineage_df.to_csv(snakemake.output.lineage_df)

# aggregate both dataframes by summing up repeating rows for variants
variants_df = variants_df.groupby(["Signatures"]).sum().reset_index()
lineage_df = lineage_df.groupby(["Signatures"]).agg("max").reset_index()

# merge variants dataframe and lineage dataframe
variants_df = variants_df.merge(lineage_df, left_on="Signatures", right_on="Signatures")

# add feature column for sorting
variants_df["Features"] = variants_df["Signatures"].str.extract(r"(.+)[:].+|\*")

# position of variant for sorting and change type
variants_df["Position"] = variants_df["Signatures"].str.extract(
    r"([0-9]+)([A-Z]+|\*)$"
)[0]
variants_df = variants_df.astype({"Position": "int64"})

# generate sorting list with correct range of features
sorterIndex = dict(zip(sorter, range(len(sorter))))
variants_df["Features_Rank"] = variants_df["Features"].map(sorterIndex)

# replace zeros with empty string and sort final DF
variants_df.replace([0, 0.0], "", inplace=True)
variants_df.sort_values(by=["Features_Rank", "Position", "VAF"], inplace=True)
variants_df.to_csv(snakemake.output[0], index=False, sep=",")
