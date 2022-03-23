import sys
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

times["s"] = pd.to_datetime(times["s"], unit="s")
# times['h:m:s'] = pd.to_datetime(times["h:m:s"], format= '%H:%M:%S')

print(times)
print(times.info())


boxplot = (
    alt.Chart(times)
    .mark_boxplot()
    .encode(
        alt.X("s:T", title=None, axis=alt.Axis(labels=False, ticks=False)),
        alt.Y("Workflow:O"),
        alt.Color("Workflow"),
    )
    .facet(column=alt.Column("Platform", title=None))
)

line = (
    alt.Chart(times)
    .mark_line(opacity=0.2)
    .encode(
        alt.X("hoursminutes(s):T", title="Wall clock time (h:m)"),
        alt.Y("# Bases"),
        alt.Color("Workflow"),
    )
)


dot = (
    alt.Chart(times)
    .mark_circle(opacity=0.4)
    .encode(
        alt.X("hoursminutes(s):T"),
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


(boxplot & dot_line).configure_legend(labelLimit=0).configure_axis(labelLimit=0).save(
    snakemake.output[0]
)

#
# .properties(
#     width=340
# )
