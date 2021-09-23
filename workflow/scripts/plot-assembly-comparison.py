sys.stderr = open(snakemake.log[0], "w")

import sys

import altair as alt
import pandas as pd
import pysam

data = pd.DataFrame()


def register_lengths(sample, file_list, state, amplicon_state, data):
    for file, assembler in zip(file_list, snakemake.params.assembler):
        if state in ("initial", "scaffolded"):
            with pysam.FastxFile(file) as infile:
                data = data.append(
                    {
                        "Sample": sample,
                        "Assembler": assembler,
                        "Amplicon": amplicon_state,
                        "length (bp)": max(len(contig.sequence) for contig in infile),
                        "State": state,
                    },
                    ignore_index=True,
                )
        else:
            quastDf = pd.read_csv(file, sep="\t")
            data = data.append(
                {
                    "Sample": sample,
                    "Assembler": assembler,
                    "Amplicon": amplicon_state,
                    "length (bp)": quastDf.loc[0, "N50"],
                    "State": "N50",
                    "Genome fraction (%)": quastDf.loc[0, "Genome fraction (%)"]
                    if "Genome fraction (%)" in quastDf.columns
                    else float("nan"),
                },
                ignore_index=True,
            )
    return data


for sample, amplicon_state in zip(
    snakemake.params.samples, snakemake.params.amplicon_state
):

    def load(inputfile, label=None):
        return register_lengths(
            sample,
            [x for x in snakemake.input.get(inputfile) if sample in x],
            label or inputfile,
            amplicon_state,
            data,
        )

    data = load("initial")
    data = load("final", "scaffolded")
    data = load("quast", "N50")

data.replace({"Amplicon": {0: "a) Shotgun", 1: "b) Amplicon"}}, inplace=True)
data.replace(
    {
        "Assembler": {
            "megahit-std": "MEGAHIT",
            "megahit-meta-large": "MEGAHIT (meta-large)",
            "megahit-meta-sensitive": "MEGAHIT (meta-sensitive)",
            "trinity": "Trinty",
            "velvet": "Velvet",
            "metaspades": "SPAdes (metaSPAdes)",
            "spades": "SPAdes",
            "coronaspades": "SPAdes (coronaSPAdes)",
            "rnaviralspades": "SPAdes (RNAviralSPAdes)",
        }
    },
    inplace=True,
)

data.to_csv(snakemake.output[1])

height = 200
width = 65

plot_bp = (
    alt.Chart(data)
    .encode(
        y=alt.Y(
            "length (bp)",
            scale=alt.Scale(domain=[0, 35000], clamp=True),
            axis=alt.Axis(tickCount=5),
        ),
        x=alt.X("State", sort=["N50", "initial", "scaffolded"]),
        color=alt.Color("Assembler", scale=alt.Scale(scheme="turbo"), legend=None),
    )
    .properties(height=height, width=width)
)

combined_bp = plot_bp.mark_point(opacity=1, filled=True) + plot_bp.mark_boxplot(
    opacity=0.8,
    box={"stroke": "black", "strokeWidth": 1, "fill": "none"},
    median={"stroke": "black", "strokeWidth": 1},
    outliers=False,
)
combined_bp = (
    combined_bp.facet(
        column=alt.Column(
            "Assembler:N",
            title="",
            header=alt.Header(labelAngle=-45, labelOrient="bottom", labelPadding=-5),
        ),
        row=alt.Row(
            "Amplicon:N",
            title="",
            sort="ascending",
            header=alt.Header(
                labelAngle=-360,
                labelAlign="center",
                labelAnchor="end",
                labelPadding=-90,
                labelOrient="left",
                labelFontWeight="bold",
                labelFontSize=12,
            ),
        ),
    )
    .configure_axisX(labelAngle=-45)
    .configure_axis(labelFontSize=12, titleFontSize=12)
    .configure_view(stroke=None)
)

plot_genome_frac = (
    alt.Chart(data)
    .mark_point(opacity=0.5)
    .encode(
        x=alt.X(
            "jitter:Q",
            stack="zero",
            title="",
            axis=alt.Axis(ticks=False, grid=False, labels=False),
            scale=alt.Scale(),
        ),
        y=alt.Y(
            "Genome fraction (%)",
            stack=None,
            title="Genome fraction (%)",
            axis=alt.Axis(tickCount=5),
        ),
        color=alt.Color("Assembler", scale=alt.Scale(scheme="turbo"), legend=None),
        column=alt.Column(
            "Assembler:N",
            title="",
            header=alt.Header(labelAngle=-45, labelOrient="bottom", labelPadding=-5),
        ),
        row=alt.Row(
            "Amplicon:N",
            title="",
            sort="ascending",
            header=alt.Header(
                labelAngle=-360,
                labelAlign="center",
                labelAnchor="end",
                labelPadding=-90,
                labelOrient="left",
                labelFontWeight="bold",
                labelFontSize=12,
            ),
        ),
    )
    .configure_axisY(grid=True)
    .configure_axis(labelFontSize=12, titleFontSize=12)
    .configure_view(stroke=None)
    .transform_calculate(jitter="20 * sqrt(-2*log(random()))*cos(2*PI*random())")
    .properties(height=height, width=50)
)


combined_bp.save(snakemake.output[0])
plot_genome_frac.save(snakemake.output[2])
