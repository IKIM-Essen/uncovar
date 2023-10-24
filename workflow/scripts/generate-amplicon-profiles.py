 # Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import pandas as pd
import numpy as np
from pysam import VariantFile
from natsort import index_natsorted, ns

bcf_in = VariantFile("/local/work/uncovar-wastewater/resources/lineage-candidate-variants/all.sorted.bcf")  # auto-detect input format

mutations = pd.DataFrame()
lineages = pd.DataFrame()

for rec in bcf_in.fetch():
    lineage_dict = {lineage: "x" for lineage in rec.info["LINEAGES"]}
    mutations = pd.concat(
        [
            mutations,
            pd.DataFrame(
                {
                    "mutation": rec.info["SIGNATURES"],
                    "position": rec.pos,
                },
        
            ),
        ],
        ignore_index=True
    )
    for signature in rec.info["SIGNATURES"]:
        lineages = pd.concat(
            [
                lineages,
                pd.DataFrame(
                    {
                        "mutation": signature,
                        **{
                            lineage.replace(".", " "): "x"
                            for lineage in rec.info["LINEAGES"]
                        },
                    },
                    index=[0],
                ),
            ],
            ignore_index=True,
        )


lineages = lineages.fillna(0)
lineages = lineages.replace({"x": 1})
lineages = (
    lineages.groupby(["mutation"])
    .agg(func={column: np.max for column in lineages.columns})
    .reset_index(drop=True)
)
lineages = lineages.replace({1: "x", 0: np.NaN})
lineages["count"] = lineages.loc[:, lineages.columns!='mutation'].count(axis=1)
print(lineages)



mutations = mutations.groupby(["mutation"]).agg(func={column: np.max for column in mutations.columns}).reset_index(drop=True)
mutations.sort_values(by=["position"], inplace=True)
mutations.to_csv("test.csv", index=False)
mutations["amplicon"] = ""
print(mutations)

mutations.set_index("mutation", inplace=True)
lineages.set_index("mutation", inplace=True)

mutations = mutations.merge(lineages, left_index=True, right_index=True)

with open("/local/work/uncovar-wastewater/D01.amplicon_stats.txt", "r") as statfile:
    for line in statfile.read().splitlines():
        if line.startswith("FREADS"):
            coverage = line.split("\t")[2:]

print(len(coverage))

amplicons = pd.read_csv("/local/work/uncovar-wastewater/resources/primer.bedpe", delimiter="\t", names=["ref_name_1", "p1_start", "p1_end", "ref_name_2", "p2_start", "p2_end"])
amplicons.drop(columns=["ref_name_1", "ref_name_2"], inplace=True)
amplicons["amp_start"] = amplicons["p1_end"] + 1
amplicons["amp_end"] = amplicons["p2_end"] - 1
amplicons["amp_len"] = amplicons["amp_end"] - amplicons["amp_start"]
amplicons["amp_cov"] = coverage

print(amplicons)

for index1, row1 in mutations.iterrows():
    for index2, row2 in amplicons.iterrows():
        if row1["position"] in range(int(row2['amp_start']), int(row2['amp_end'])):
            if mutations.loc[index1]["amplicon"] == "":
                mutations.at[index1, "amplicon"] = str(index2)
            else:
                mutations.at[index1, "amplicon"] += "," + str(index2)

all_columns = mutations.columns
first_columns = ["position", "amplicon", "count"]
rest_columns = [item for item in all_columns.sort_values() if item not in first_columns]
mutations = mutations[[*first_columns, *rest_columns]]


mutations = (mutations.drop('amplicon', axis=1)
            .join
            (
            mutations.amplicon
            .str
            .split(",", expand=True)
            .stack()
            .reset_index(drop=True, level=1)
            .rename('amplicon')           
            ))

mutations = mutations[[*first_columns, *rest_columns]]
mutations.sort_values(by=["position", "amplicon"], inplace=True)
mutations.to_csv("amplicons.tsv", index=True, sep="\t")



amplicon_profiles = mutations


# amplicon_profiles.set_index("mutation", inplace=True)
for index, row in amplicon_profiles.iterrows():
    for column in rest_columns:
        if row[column] == "x":
            amplicon_profiles.at[index, column] = index

amplicon_profiles.reset_index(drop=True, inplace=True)
amplicon_profiles.drop(columns=["position", "count"], inplace=True)
amplicon_profiles = amplicon_profiles.replace({np.NaN : ""})
amplicon_profiles = amplicon_profiles.groupby(["amplicon"]).agg(lambda x: ' '.join(set(x)))
amplicon_profiles.sort_index(key= lambda x: np.argsort(index_natsorted(amplicon_profiles.index)), inplace=True)
print(amplicon_profiles)

amplicon_profiles.to_csv("amp_profiles.csv", index=True)
# pos = 21764
# test = amplicons[(pos >= amplicons['amp_start']) & (pos <= amplicons['amp_end'])]
# print(test)

# print(mutations)

# import sys

# sys.stderr = open(snakemake.log[0], "w")

# import altair as alt
# import pandas as pd
# import pysam
# from intervaltree import IntervalTree

# # read primer bedpe to df
# PRIMER = pd.read_csv(snakemake.params.get("bedpe", ""), delimiter="\t", header=None)
# PRIMER.drop(PRIMER.columns[[0, 3]], axis=1, inplace=True)
# PRIMER.columns = ["p1_start", "p1_end", "p2_start", "p2_end"]

# # convert df to interval trees
# primer_intervals = IntervalTree()
# no_primer_intervals = IntervalTree()
# for index, row in PRIMER.iterrows():
#     primer_intervals[row["p1_start"] : row["p2_end"] + 1] = (
#         row["p1_start"],
#         row["p2_end"] + 1,
#     )
#     no_primer_intervals[row["p1_end"] + 1 : row["p2_start"]] = (
#         row["p1_end"] + 1,
#         row["p2_start"],
#     )


# def iter_with_samples(inputfiles):
#     return zip(snakemake.params.samples, inputfiles)


# def count_intervals(file):
#     with pysam.AlignmentFile(file, "rb") as bam:
#         counter_primer = 0
#         counter_no_primer = 0
#         counter_primer_within = 0
#         counter_no_primer_within = 0
#         counter_nothing = 0
#         mate_pair_intervals = {}
#         for read in bam.fetch():
#             if not mate_pair_intervals.get(read.query_name):
#                 mate_pair_intervals[read.query_name] = [read.reference_start]
#             else:
#                 mate_pair_intervals[read.query_name].append(read.reference_end)
#         for pair in mate_pair_intervals:
#             if (
#                 len(mate_pair_intervals[pair]) > 1
#                 and mate_pair_intervals[pair][0] != None
#                 and mate_pair_intervals[pair][1] != None
#             ):
#                 if primer_intervals.envelop(
#                     mate_pair_intervals[pair][0], mate_pair_intervals[pair][1] + 1
#                 ):
#                     if (
#                         sorted(
#                             primer_intervals.envelop(
#                                 mate_pair_intervals[pair][0],
#                                 mate_pair_intervals[pair][1] + 1,
#                             )
#                         )[0].begin
#                         == mate_pair_intervals[pair][0]
#                         and sorted(
#                             primer_intervals.envelop(
#                                 mate_pair_intervals[pair][0],
#                                 mate_pair_intervals[pair][1] + 1,
#                             )
#                         )[0].end
#                         == mate_pair_intervals[pair][1] + 1
#                     ):
#                         counter_primer += 1
#                     else:
#                         counter_primer_within += 1
#                 elif no_primer_intervals.envelop(
#                     mate_pair_intervals[pair][0] + 1, mate_pair_intervals[pair][1]
#                 ):
#                     if (
#                         sorted(
#                             no_primer_intervals.envelop(
#                                 mate_pair_intervals[pair][0] + 1,
#                                 mate_pair_intervals[pair][1],
#                             )
#                         )[0].begin
#                         == mate_pair_intervals[pair][0] + 1
#                         and sorted(
#                             no_primer_intervals.envelop(
#                                 mate_pair_intervals[pair][0] + 1,
#                                 mate_pair_intervals[pair][1],
#                             )
#                         )[0].end
#                         == mate_pair_intervals[pair][1]
#                     ):
#                         counter_no_primer += 1
#                     else:
#                         counter_no_primer_within += 1
#                 else:
#                     counter_nothing += 1
#             else:
#                 counter_nothing += 1
#         counters = pd.DataFrame(
#             {
#                 "n_count": [
#                     counter_primer,
#                     counter_primer_within,
#                     counter_no_primer,
#                     counter_no_primer_within,
#                     counter_nothing,
#                 ],
#                 "class": [
#                     "uncut primer exact",
#                     "uncut primer within",
#                     "cut primer exact",
#                     "cut primer within",
#                     "no mathing win",
#                 ],
#             }
#         )
#         return counters


# def plot_classes(counters):
#     bars = (
#         alt.Chart(counters)
#         .mark_bar()
#         .encode(
#             y="class",
#             x="n_count",
#             row=alt.Row("sample:N"),
#             column=alt.Column("state:N", sort="descending"),
#         )
#     )
#     text = bars.mark_text(
#         align="left",
#         baseline="middle",
#         dx=3,  # Nudges text to right so it doesn't appear on top of the bar
#     ).encode(text="n_count", row=alt.Row("sample:N"), column=alt.Column("state:N"))
#     return bars, text


# all_df = pd.DataFrame()
# for sample, file in iter_with_samples(snakemake.input.unclipped):
#     counts_before = count_intervals(file)
#     counts_before["sample"] = sample
#     counts_before["state"] = "before"
#     all_df = pd.concat([all_df, counts_before], ignore_index=True)

# for sample, file in iter_with_samples(snakemake.input.clipped):
#     counts_after = count_intervals(file)
#     counts_after["sample"] = sample
#     counts_after["state"] = "after"
#     all_df = pd.concat([all_df, counts_after], ignore_index=True)

# bars, text = plot_classes(all_df)

# (bars).properties(title="Amplicon matching").save(snakemake.output.plot)
