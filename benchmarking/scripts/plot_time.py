import sys
from cgi import print_environ
from time import time
from turtle import width

sys.stderr = open(snakemake.log[0], "w")

import altair as alt
import pandas as pd

WORKFLOWS = {
    "ncov2019-artic-nf-": "ncov2019",
    "ncov2019-artic-nf-medaka": "ncov2019 (medaka)",
    "ncov2019-artic-nf-nanopolish": "ncov2019 (nanopolish)",
    "artic-medaka": "ARTIC (medaka)",
    "artic-nanopolish": "ARTIC (nanopolish)",
    "nf-core-viralrecon-medaka": "viralrecon (medaka)",
    "nf-core-viralrecon-nanopolish": "viralrecon (nanopolish)",
    "uncovar": "UnCoVar",
    "havoc": "HAVoC",
    "covpipe": "CoVpipe",
    "snakelines": "SnakeLines",
    "v-pipe": "V-pipe",
    "signal": "SIGNAL",
}

PLATTFORM = {"illumina": "Illumina Workflows", "nanopore": "Nanopore Workflows"}

times = pd.read_csv(snakemake.input.times, sep="\t")

size = []
for path in snakemake.input.sizes:
    size.append(pd.read_csv(path, sep="\t"))

size = pd.concat(size, axis=0, ignore_index=True)

times = times.merge(size, how="left")

times["Workflow"] = times["Workflow"].replace(WORKFLOWS)
times["Platform"] = times["Platform"].replace(PLATTFORM)

times["s"] = pd.to_timedelta(times["s"], unit="s")
times["s"] = times["s"].dt.total_seconds().div(60)
times["cpu_time"] = pd.to_timedelta(times["cpu_time"], unit="s")
times["cpu_time"] = times["cpu_time"].dt.total_seconds().div(60)


def boxplot(x):
    return (
        alt.Chart(times)
        .mark_boxplot()
        .encode(
            alt.X(x, title=None, axis=alt.Axis(labels=False, ticks=False)),
            alt.Y("Workflow:O"),
            alt.Color("Workflow"),
        )
        .facet(column=alt.Column("Platform", title=None))
    )


def dot_line(x, title):
    line = (
        alt.Chart(times)
        .mark_line(opacity=0.2)
        .encode(
            alt.X(x, title=title),
            alt.Y("# Bases"),
            alt.Color("Workflow"),
        )
    )

    dot = (
        alt.Chart(times)
        .mark_circle(opacity=0.4)
        .encode(
            alt.X(x),
            alt.Y("# Bases"),
            alt.Color(
                "Workflow",
                legend=alt.Legend(
                    orient="top", direction="vertical", columns=4, symbolLimit=0
                ),
                title=None,
            ),
        )
    )

    dot_line = (line + dot).facet(
        column=alt.Column("Platform", title=None, header=alt.Header(labels=False))
    )

    return dot_line


(
    # boxplot("s:Q") &
    # dot_line("s:Q",  "Wall clock time (min)")
    boxplot("cpu_time:Q")
    & dot_line("cpu_time:Q", "CPU time (min)")
).configure_legend(labelLimit=0).configure_axis(labelLimit=0).save(snakemake.output[0])
