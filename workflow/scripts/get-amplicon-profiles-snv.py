# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import pandas as pd
import numpy as np
from pysam import VariantFile
from natsort import index_natsorted, ns
from more_itertools import powerset

bcf_in = VariantFile("/local/work/uncovar-wastewater/resources/lineage-candidate-variants/all.sorted.bcf")  # auto-detect input format

mutations = pd.DataFrame()
lineages = pd.DataFrame()

for rec in bcf_in.fetch():
    for signature in rec.info["SIGNATURES"]:
        for lineage in rec.info["LINEAGES"]:
            lineages = pd.concat(
                [
                    lineages,
                    pd.DataFrame(
                        {
                            "mutation": signature,
                            "lineage": lineage.replace(".", " "),
                            "position": rec.pos,
                        },
                        index=[0],
                    ),
                ],
                ignore_index=True,
            )

lineages.drop_duplicates(inplace=True)
lineages.reset_index(inplace=True, drop=True)

voc_names = pd.read_csv("/local/work/uncovar-wastewater/resources/variants-rename.csv", header=0)
print(voc_names)
newer_omicrons = voc_names[
    (voc_names["WHO Label"] == "Omicron") & 
    (voc_names["Nextstrain Clade"] != "21K Omicron") & 
    (voc_names["Nextstrain Clade"] != "21L Omicron")
    ]["Nextstrain Clade"].to_list()

exclude_list = ["20A S 126A", "21J Delta", "21I Delta"]
# lineages["lineage"].isin(newer_omicrons).index.drop()
# lineages.drop(lineages[lineages["lineage"].isin(newer_omicrons)].index, inplace=True)
lineages.drop(lineages[lineages["lineage"].isin(exclude_list)].index, inplace=True)
# lineages[lineages.groupby('mutation').mutation.transform(len) == 1]
print(lineages[lineages.groupby('mutation').mutation.transform(len) == 1])
print(len(lineages.mutation.unique()))
# lineages.to_csv("snv.csv")

# files = [
#     "/local/work/uncovar-wastewater/results/2022-06-30/lineage-variant-report/sewerA-04.csv",
#     "/local/work/uncovar-wastewater/results/2022-06-30/lineage-variant-report/sewerA-05.csv",
#     "/local/work/uncovar-wastewater/results/2022-06-30/lineage-variant-report/sewerA-06.csv",
#     "/local/work/uncovar-wastewater/results/2022-06-30/lineage-variant-report/sewerB-04.csv",
#     "/local/work/uncovar-wastewater/results/2022-06-30/lineage-variant-report/sewerB-05.csv",
#     "/local/work/uncovar-wastewater/results/2022-06-30/lineage-variant-report/sewerB-06.csv",
# ]

# # # samples = ["D01"]
# samples = [
#     "sewerA-04",
#     "sewerA-05",
#     "sewerA-06",
#     "sewerB-04",
#     "sewerB-05",
#     "sewerB-06",
# ]

data = pd.DataFrame()
# data = pd.read_csv(snakemake.input.var_df, usecols=["Mutations","Probability","Frequency","ReadDepth"])

vocs_included = {}
# for file, sample in zip(files, samples):
for file, sample in zip(snakemake.input.csv, snakemake.params.sample):
    data = pd.concat(
        [
            data,
            pd.read_csv(file, usecols=["Mutations","Probability","Frequency","ReadDepth"])
        ]
    )

    data.set_index("Mutations", inplace=True)
    data.index.rename("mutation", inplace=True)
    # mutations.
    data.drop(["Lineage", "Similarity"], inplace=True)
    data.drop(["ORF1a:L3674-", "ORF1a:S2083-"], inplace=True)

    filter_data = data[(data["Probability"] > 0.95) & (data["Frequency"] * data["ReadDepth"] > 10)].copy(deep=True)

    # print(filter_data)

    rename_dict = voc_names.set_index('Nextstrain Clade').to_dict()['Pango Lineage']
    # print(rename_dict)
    # lineages["lineage"].replace(rename_dict, inplace=True)

    out_list = []
    for mutation in filter_data.index:
        eindeut_lineage = lineages[lineages["mutation"] == mutation]["lineage"].to_list()
        eindeut_lineage = [x for x in eindeut_lineage if x != "nan"]
        if len(eindeut_lineage) <= 1:
            out_list.append("/".join(eindeut_lineage))
        print(mutation, lineages[lineages["mutation"] == mutation]["lineage"].to_list())
    out_list = list(dict.fromkeys(out_list))
    out_list = [rename_dict[x] if  x in rename_dict else x for x in out_list]
    out_list = [x for x in out_list if x != ""]
    site, week = sample.split("-")
    # site = sample[:-2]
    # week = sample[-2:]
    if week not in vocs_included:
        vocs_included[week] = {}
    vocs_included[week][site] = out_list

print(vocs_included)

vocs = pd.DataFrame()
for week in vocs_included:
    for site in vocs_included[week]:
        vocs.at[week, site] = " ".join(vocs_included[week][site])

print(vocs)
vocs.to_csv(snakemake.output.ampprofile)