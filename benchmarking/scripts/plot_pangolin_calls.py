import sys
from time import altzone

sys.stderr = open(snakemake.log[0], "w")

import altair as alt
import pandas as pd

source = pd.read_csv(snakemake.input[0])

rectangles = (
    alt.Chart()
    .mark_rect()
    .encode(
        alt.X("Workflow:N"),
        alt.Y("Sample:N"),
        alt.Color("lineage:N"),
    )
)

text = (
    alt.Chart()
    .mark_text()
    .encode(alt.X("Workflow:N"), alt.Y("Sample:N"), alt.Text("lineage:N"))
)

(rectangles + text).facet(column="Platform:N", data=source).save(snakemake.output[0])
