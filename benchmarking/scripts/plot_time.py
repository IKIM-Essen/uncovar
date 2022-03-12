import sys

sys.stderr = open(snakemake.log[0], "w")

import altair as alt
import pandas as pd

times = pd.read_csv(snakemake.input.times, sep="\t")

size = []
for path in snakemake.input.sizes:
    size.append(pd.read_csv(path, sep="\t"))

size = pd.concat(size, axis=0, ignore_index=True)

times = times.merge(size, how="left")

print(times)
print(times.info())

alt.Chart(times).mark_point().encode(
    alt.X("minutes(s)", title="Minutes"), alt.Y("# Bases"), alt.Color("Workflow")
).save(snakemake.output[0])
