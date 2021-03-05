sys.stderr = open(snakemake.log[0], "w")
min_coverage = snakemake.params.get("min_coverage", "")

import pandas as pd
import altair as alt


def plot_coverage(sm_input, sm_output, min_coverage):
    coverage = pd.read_csv(sm_input, sep="\t")

    coverage.rename(
        columns={coverage.columns[2]: "Coverage", "POS": "Pos"}, inplace=True
    )

    coverage["Sample"] = coverage["#CHROM"].apply(lambda x: str(x).split(".")[0])
    coverage["# Coverage"] = coverage.Coverage.apply(
        lambda x: f"< {min_coverage}"
        if int(x) < int(min_coverage)
        else f">= {min_coverage}"
    )

    # max_x_pos = 30000 if coverage.Pos.max() < 30000 else coverage.Pos.max()
    max_x_pos = coverage.Pos.max()

    if len(coverage) > 0:
        alt.Chart(coverage).mark_bar().encode(
            x=alt.X("Pos:Q", scale=alt.Scale(domain=(0, max_x_pos), nice=False)),
            y=alt.Y(
                "Coverage", scale=alt.Scale(domain=[0, coverage.Coverage.max() * 1.05])
            ),
            column=alt.Column("Sample:N"),
            color=alt.Color(
                "# Coverage",
                scale=alt.Scale(
                    domain=[f"< {min_coverage}", f">= {min_coverage}"],
                    range=["indianred", "lightgreen"],
                ),
            ),
        ).properties(width=600, height=150).save(sm_output)
    else:
        alt.Chart(coverage).mark_bar().encode().properties(width=600, height=150).save(
            sm_output
        )


plot_coverage(snakemake.input[0], snakemake.output[0], min_coverage)
