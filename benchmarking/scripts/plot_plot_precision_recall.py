import sys

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

value_vars = ["METRIC-Recall", "METRIC-Precision"]
id_vars = metrics.columns.tolist()

id_vars.remove("Unnamed: 0")

for value_var in value_vars:
    id_vars.remove(value_var)

metrics = metrics.melt(
    id_vars=id_vars, value_vars=value_vars, var_name="Metric", value_name="Value"
)

# metrics.columns = [col.removeprefix("METRIC-") for col in metrics.columns]
metrics["Metric"] = metrics["Metric"].str.removeprefix("METRIC-")

for key, value in WORKFLOWS.items():
    metrics["Workflow"] = metrics["Workflow"].str.replace(key, value)

for key, value in PLATTFORM.items():
    metrics["Platform"] = metrics["Platform"].str.replace(key, value)

metrics["Workflow"] = metrics["Workflow"] + " (" + metrics["Mode"].fillna("") + ")"
metrics["Workflow"] = metrics["Workflow"].str.removesuffix(" ()")

bars = (
    alt.Chart(metrics)
    .mark_bar()
    .encode(
        alt.X("mean(Value):Q"), alt.Y("Workflow:N", title=None), alt.Color("Type:N")
    )
)

error_bars = (
    alt.Chart()
    .mark_errorbar(extent="ci")
    .encode(
        alt.X("Value:Q"),
        alt.Y("Workflow:N"),
    )
)

plot = alt.layer(bars, error_bars, data=metrics).facet(
    column="Platform:N", row="Metric:N"
)


plot.save(snakemake.output[0])
