# TODO check whether this is still in use

import sys

sys.stderr = open(snakemake.log[0], "w")

MIXTURE_PREFIX = snakemake.params.get("prefix", "")
MIXTURE_PART_INDICATOR = snakemake.params.get("separator", "")
MIXTURE_PERCENTAGE_INDICATOR = snakemake.params.get("percentage", "")

import pandas as pd
import altair as alt
import re


def regex_per_line(content, pattern):
    match = re.search(pattern, content)
    if match:
        return int(match.group("share")) / 100
    else:
        return 0


def plot_dependency_of_pangolin_call(sm_input, sm_output):
    # aggregate pangolin outputs
    all_sampes = pd.DataFrame()
    for input in sm_input:
        pangolin_output = pd.read_csv(input)
        pangolin_output["mixture_content"] = input.split(MIXTURE_PREFIX, 1)[-1].split(
            "."
        )[0]
        all_sampes = all_sampes.append(pangolin_output, ignore_index=True)

    all_sampes["mixture_content"] = all_sampes["mixture_content"].str.replace("-", ".")

    # get share of called lineage
    all_sampes["regex"] = (
        MIXTURE_PART_INDICATOR
        + all_sampes["lineage"]
        + MIXTURE_PERCENTAGE_INDICATOR
        + "(?P<share>\d.)"
    )
    all_sampes["share_of_called_lineage"] = all_sampes[
        ["mixture_content", "regex"]
    ].apply(lambda x: regex_per_line(*x), axis=1)
    all_sampes.drop(columns=["regex"], inplace=True)
    all_sampes["correct_call"] = all_sampes["share_of_called_lineage"].apply(
        lambda x: "Yes" if x > 0 else "No"
    )

    # plot histogram
    source = all_sampes.copy()
    source.rename(
        columns={
            "share_of_called_lineage": "Share of called lineage in sample",
            "correct_call": "Lineage in sample",
        },
        inplace=True,
    )

    histogram = (
        alt.Chart(source)
        .mark_bar()
        .encode(
            alt.X(
                "Share of called lineage in sample:Q",
                bin=alt.Bin(extent=[0, 1], step=0.1),
                axis=alt.Axis(format="%"),
            ),
            y="count()",
            color=alt.Color(
                "Lineage in sample",
                scale=alt.Scale(scheme="tableau20"),
                sort="descending",
            ),
        )
        .configure_legend(orient="top")
    )

    histogram.save(sm_output)


if __name__ == "__main__":
    plot_dependency_of_pangolin_call(snakemake.input, snakemake.output[0])
