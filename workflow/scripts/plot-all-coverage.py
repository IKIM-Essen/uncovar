import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import altair as alt


def plot_coverage(sm_input, sm_output, min_coverage):

    coverage = pd.DataFrame()
    for sample in sm_input:
        sample_df = pd.read_csv(sample, sep="\t")

        sample_df.rename(
            columns={sample_df.columns[2]: "Coverage", "POS": "Pos"}, inplace=True
        )

        sample_df["Sample"] = sample_df["#CHROM"].apply(lambda x: str(x).split(".")[0])

        coverage = coverage.append(sample_df, ignore_index=True)

    coverage["# Coverage"] = coverage.Coverage.apply(
        lambda x: f"< {min_coverage}"
        if int(x) < int(min_coverage)
        else f">= {min_coverage}"
    )

    max_y_pos = 100
    max_x_pos = coverage.Pos.max()

    coverage["Coverage"] = coverage["Coverage"].apply(
        lambda x: max_y_pos if x > max_y_pos else x
    )

    if len(coverage) > 0:
        alt.Chart(coverage).mark_bar().encode(
            x=alt.X("Pos:Q", scale=alt.Scale(domain=(0, max_x_pos), nice=False)),
            y=alt.Y("Coverage", scale=alt.Scale(domain=[0, max_y_pos])),
            row=alt.Row("Sample:N"),
            color=alt.Color(
                "# Coverage",
                scale=alt.Scale(
                    domain=[f"< {min_coverage}", f">= {min_coverage}"],
                    range=["indianred", "lightgreen"],
                ),
            ),
        ).properties(width=1200, height=150).interactive().save(sm_output)
    else:
        alt.Chart(coverage).mark_bar().encode().properties(
            width="container", height=150
        ).save(sm_output)


plot_coverage(snakemake.input, snakemake.output[0], snakemake.params.min_coverage)
