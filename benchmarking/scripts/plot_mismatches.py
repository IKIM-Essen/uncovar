import sys
from os import sep
from turtle import color

import altair as alt
import pandas as pd

# sys.stderr = open(snakemake.log[0], "w")


SUFFIXES = ["-medaka", "-nanopolish"]

WORKFLOWS = {
    "ncov2019-artic-nf": "ncov2019",
    "artic": "ARTIC",
    "nf-core-viralrecon": "viralrecon",
    "uncovar": "UnCoVar",
    "havoc": "HAVoC",
    "covpipe": "CoVpipe",
    "snakelines": "SnakeLines",
    "v-pipe": "V-pipe",
    "signal": "SIGNAL",
}

PLATTFORM = {"illumina": "Illumina", "nanopore": "Nanopore"}

metrics = pd.read_csv(snakemake.input[0], sep="\t")

metrics.columns = [col.replace(".", "-") for col in metrics.columns]

for key, value in WORKFLOWS.items():
    metrics["Workflow"] = metrics["Workflow"].str.replace(key, value)

for key, value in PLATTFORM.items():
    metrics["Platform"] = metrics["Platform"].str.replace(key, value)

metrics["Workflow"] = metrics["Workflow"] + " (" + metrics["Mode"].fillna("") + ")"
metrics["Workflow"] = metrics["Workflow"].str.removesuffix(" ()")

metrics.drop(
    columns=[
        "METRIC-Recall",
        "METRIC-Precision",
        "METRIC-Frac_NA",
        "METRIC-F1_Score",
        "TRUTH-TOTAL-TiTv_ratio",
        "QUERY-TOTAL-TiTv_ratio",
        "TRUTH-TOTAL-het_hom_ratio",
        "QUERY-TOTAL-het_hom_ratio",
        "Unnamed: 0",
    ],
    inplace=True,
)

metrics = metrics.groupby(by=["Workflow", "Platform", "Type"]).sum()
metrics["Recall"] = metrics["TRUTH-TP"] / (metrics["TRUTH-TP"] + metrics["TRUTH-FN"])
metrics["Precision"] = metrics["TRUTH-TP"] / (metrics["TRUTH-TP"] + metrics["QUERY-FP"])
metrics.reset_index(inplace=True)

print(metrics)

value_vars = ["Recall", "Precision"]

id_vars = metrics.columns.tolist()
for value_var in value_vars:
    id_vars.remove(value_var)

metrics = metrics.melt(id_vars=id_vars, value_vars=value_vars, var_name="Metric")

metrics = metrics.drop(
    columns=[
        "TRUTH-TOTAL",
        "TRUTH-TP",
        "TRUTH-FN",
        "QUERY-TOTAL",
        "QUERY-FP",
        "Metric",
        "value",
        "QUERY-UNK",
    ]
).melt(
    id_vars=["Workflow", "Platform", "Type"],
    value_vars=["FP-gt", "FP-al"],
    var_name="Mismatches",
    value_name="Number",
)

metrics["Mismatches"].replace(
    {
        "FP-gt": "Genotype",
        "FP-al": "Allelic",
    },
    inplace=True,
)


def plot_mismatches(data, platform):
    return (
        alt.Chart()
        .mark_bar()
        .encode(
            alt.Y("Number:Q"),
            alt.X("Mismatches:N", title=None),
            alt.Color("Mismatches:N"),
        )
        + alt.Chart()
        .mark_text(
            dy=-5,
        )
        .encode(
            alt.Y("Number:Q"), alt.X("Mismatches:N", title=None), alt.Text("Number:Q")
        )
    ).facet(
        column=alt.Column("Workflow:N", title=f"Variant Calls on {platform} Workflows"),
        row=alt.Row("Type:N", title=None),
        data=data,
    )


ill = metrics.loc[metrics["Platform"] == "Illumina"]
ont = metrics.loc[metrics["Platform"] == "Nanopore"]

metrics.to_csv(snakemake.output.data, sep="\t", index=False)

alt.vconcat(plot_mismatches(ill, "Illumina"), plot_mismatches(ont, "Nanopore")).save(
    snakemake.output.plot
)
