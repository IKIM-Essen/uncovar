# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys
import re

sys.stderr = open(snakemake.log[0], "w")
# sys.stdout = open(snakemake.log[0], "a")

import pandas as pd
import pysam
import gffutils

def phred_to_prob(phred):
    if phred is None:
        return 0
    return 10 ** (-phred / 10)

def has_numbers(inputString):
    return any(char.isdigit() for char in inputString)

variants = pd.DataFrame()
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
            variants = variants.append(

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
                ignore_index=True
            )
            lineage_df = lineage_df.append(
                {
                    "Signatures": signatures,
                    **lineage_dict
                },
                ignore_index=True
            )


sigs = list(variants["Signatures"])
print(sigs)
print(len(sigs))
sigs = list(dict.fromkeys(sigs))
print(sorted(sigs))
print(len(sigs))
variants.to_csv(snakemake.output.lineage_df)
variants = variants.groupby(["Signatures"]).sum().reset_index()
print(variants)
# print(variants)
# calculation = lambda value: 
# variants = variants.groupby(["Signatures"]).agg(sum_VAF=("VAF", "sum"), sum_PROB=("Prob_clonal", "sum")).reset_index()
# lineage_df.drop_duplicates(inplace=True)
lineage_df = lineage_df.groupby(["Signatures"]).agg("max").reset_index()
lineage_df.sort_values("Signatures", inplace=True)
variants = variants.merge(lineage_df, left_on="Signatures", right_on="Signatures")
# lineage_df.to_csv(snakemake.output.lineage_df)
sigs = list(lineage_df["Signatures"])
print(sigs)
print(len(sigs))
sigs = list(dict.fromkeys(sigs))
print(sorted(sigs))
print(len(sigs))
variants["Features"] = variants["Signatures"].str.extract(r'(.+)[:].+|\*')
# print(no_number)
variants.sort_values(by=["Signatures"], inplace=True)

# variants["is AA"] = variants["Signatures"].str.contains(":")
variants["Position"] = variants["Signatures"].str.extract(r'([0-9]+)([A-Z]+|\*)$')[0]
variants = variants.astype({'Position':'int64'})

sorterIndex = dict(zip(sorter, range(len(sorter))))
print(sorterIndex)
variants["Features_Rank"] = variants["Features"].map(sorterIndex)

print(variants.dtypes)
variants.replace([0, 0.0], '', inplace=True)
# variants.sort_values(by=["Signatures_Rank", "Position", "is AA"], ascending=[True, True, False], inplace=True)
variants.sort_values(by=["Features_Rank", "Position"], inplace=True)
variants.to_csv(snakemake.output[0], index=False, sep=",")

    