import sys
from turtle import color

import altair as alt
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")


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

# SNP           SNP or MNP variants. We count single nucleotides that have changed
# INDEL         Indels and complex variants
# TRUTH.TOTAL	Total number of truth variants
# TRUTH.TP	    Number of true-positive calls in truth representation (counted via the truth sample column)
# TRUTH.FN	    Number of false-negative calls = calls in truth without matching query call
# QUERY.TOTAL	Total number of query calls
# QUERY.TP	    Number of true positive calls in query representation (counted via the query sample column)
# QUERY.FP	    Number of false-positive calls in the query file (mismatched query calls within the confident regions)
# QUERY.UNK	    Number of query calls outside the confident regions
# FP.gt	        Number of genotype mismatches (alleles match, but different zygosity)
# FP.al	        Number of allele mismatches (variants matched by position and not by haplotype)

# Recall = TP/(TP+FN)
# Precision = TP/(TP+FP)
# Frac_NA = UNK/total(query)
# F1_Score = 2 * Precision * Recall / (Precision + Recall)

metrics = metrics.groupby(by=["Workflow", "Platform", "Type"]).sum()

metrics["TP_no_gt"] = metrics["QUERY-TOTAL"] - metrics["QUERY-FP"] + metrics["FP-gt"]
metrics["T"] = metrics["TP_no_gt"] + metrics["TRUTH-FN"]
metrics["Precision"] = metrics["TP_no_gt"] / (metrics["TP_no_gt"] + metrics["TRUTH-FN"])
metrics["Recall"] = metrics["TP_no_gt"] / (metrics["TP_no_gt"] + metrics["QUERY-FP"])

metrics.reset_index(inplace=True)

value_vars = ["Recall", "Precision"]

id_vars = metrics.columns.tolist()
for value_var in value_vars:
    id_vars.remove(value_var)

metrics = metrics.melt(id_vars=id_vars, value_vars=value_vars, var_name="Metric")


def get_number(metric, value, tp, qt, t):
    if metric == "Recall":  # TP_no_gt / QUERY-TOTAL
        return f"{round(value,2)} ({int(tp)}/{int(qt)})"
    elif metric == "Precision":  # TP_no_gt / T
        return f"{round(value,2)} ({int(tp)}/{int(t)})"


metrics["Numbers"] = metrics.apply(
    lambda row: get_number(
        metric=row["Metric"],
        value=row["value"],
        tp=row["TP_no_gt"],
        qt=row["QUERY-TOTAL"],
        t=row["T"],
    ),
    axis=1,
)

metrics.to_csv(snakemake.output.data, sep="\t", index=False)


def barplot(platform):
    return (
        alt.Chart()
        .mark_bar()
        .encode(
            alt.X("value:Q", title=None),
            alt.Y("Metric:N", title=None),
            alt.Color("Metric:N", legend=None),
        )
    )


def plot_numbers(platform):
    return (
        alt.Chart()
        .mark_text(
            color="black",
            align="left",
            baseline="middle",
            dx=4,
        )
        .encode(
            alt.X("value:Q"),
            alt.Y("Metric:N", title=None),
            alt.Text("Numbers"),
        )
    )


def faceted(data, platform):
    chart = barplot(platform) + plot_numbers(platform)
    return chart.facet(
        row=alt.Row(
            "Workflow:N", title=None, header=alt.Header(labelAngle=0, labelAlign="left")
        ),
        column=alt.Column("Type:N", title=f"Variant Calls on {platform} Workflows"),
        data=data,
    )


ill = metrics.loc[metrics["Platform"] == "Illumina"]
ont = metrics.loc[metrics["Platform"] == "Nanopore"]

alt.vconcat(faceted(ill, "Illumina"), faceted(ont, "Nanopore"), data=metrics).save(
    snakemake.output.plot
)
