import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import altair as alt


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
    pangolin_calls["lineage_count"] = pangolin_calls.groupby("lineage", as_index=False)[
        "lineage"
    ].transform(lambda s: s.count())

    # mask low occurrences
    pangolin_calls.loc[
        pangolin_calls["lineage_count"] < 10, "lineage"
    ] = "other (< 10 occ.)"

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
