# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys

sys.stderr = open(snakemake.log[0], "w")

import altair as alt
import pandas as pd


def plot_lineages_over_time(sm_input, sm_output, dates, sm_output_table):
    pangolin_outputs = []
    for call, date in zip(sm_input, dates):
        pangolin_call = pd.read_csv(call)
        pangolin_call["date"] = date
        pangolin_outputs.append(pangolin_call)

    pangolin_calls = pd.concat(pangolin_outputs, axis=0, ignore_index=True)

    # write out as table
    pangolin_calls.to_csv(sm_output_table)

    pangolin_calls = pangolin_calls[pangolin_calls["lineage"] != "None"]

    # get occurrences
    if len(pangolin_calls) > 0:
        pangolin_calls["lineage_count"] = pangolin_calls.groupby(
            "lineage", as_index=False
        )["lineage"].transform(lambda s: s.count())
    else:
        pangolin_calls["lineage_count"] = pd.Series()

    # mask low occurrences
    print(pangolin_calls["lineage"].value_counts())
    df = pd.DataFrame(pangolin_calls["lineage"].value_counts())
    df.sort_values(by=["lineage"])
    if len(df.index) > 10:
        pangolin_calls.loc[
            ~pangolin_calls["lineage"].isin(df.head(10).index), "lineage"
        ] = "other occ."
    pangolin_calls.rename(columns={"lineage": "Lineage", "date": "Date"}, inplace=True)

    area_plot = (
        alt.Chart(pangolin_calls)
        .mark_bar(opacity=0.8)
        .encode(
            x=alt.X("Date:O"),
            y=alt.Y(
                "count()",
                stack="normalize",
                axis=alt.Axis(format="%"),
                title="Fraction in Run",
            ),
            stroke="Lineage",
            color=alt.Color(
                "Lineage",
                scale=alt.Scale(scheme="tableau10"),
                legend=alt.Legend(orient="top"),
            ),
        )
    ).properties(width=800)

    area_plot.save(sm_output)


plot_lineages_over_time(
    snakemake.input, snakemake.output[0], snakemake.params.dates, snakemake.output[1]
)
