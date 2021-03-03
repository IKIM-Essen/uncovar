sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import altair as alt


def plot_coverage(sm_input, sm_output):
    coverage = pd.read_csv(sm_input, sep="\t")
    coverage.rename(
        columns={coverage.columns[2]: "Coverage", "POS": "Pos"}, inplace=True
    )
    coverage["Sample"] = coverage["#CHROM"].apply(lambda x: str(x).split(".")[0])

    alt.Chart(coverage).mark_bar().encode(
        x=alt.X("Pos:Q", scale=alt.Scale(domain=(0, coverage.Pos.max()), nice=False)),
        y=alt.Y(
            "Coverage", scale=alt.Scale(domain=[0, coverage.Coverage.max() * 1.05])
        ),
        column=alt.Column("Sample:N"),
    ).properties(width=600, height=150).save(snakemake.output[0])


plot_coverage(snakemake.input[0], snakemake.output[0])
