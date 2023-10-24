 # Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import pandas as pd
import altair as alt

primer_bed = pd.read_csv(snakemake.input.bedpe, delimiter="\t", header=None)
primer_count = len(primer_bed.index)

print(primer_count)


amp_num = [*range(1, primer_count + 1)]

# names = pd.read_csv(snakemake.input.names, header=0, index_col="processing_name")
# weeks = pd.read_csv(snakemake.input.weeks, header=0, index_col="number", )

amp_cov = pd.DataFrame()
amp_cov["amplicon"] = amp_num


for sample, file in zip(snakemake.params.samples, snakemake.input.amp_stats):
    with open(file, "r") as statfile:
        # sample_rename = f"{names.loc[sample[0]]['screen_name']} {weeks.loc[int(sample[1:])]['week_number']}"
        for line in statfile.read().splitlines():
            if line.startswith("FDEPTH"):
                coverage = line.split("\t")[2:]
                amp_cov[sample] = [float(x) for x in coverage]
                # amp_cov[sample_rename] = [float(x) for x in coverage]

amp_cov.to_csv(snakemake.output.stats)
amp_cov = amp_cov.melt('amplicon', var_name='sample', value_name='coverage')
print(amp_cov["coverage"].dtypes)
amp_cov = amp_cov.round({"coverage": 0})

amp_cov.sort_values(by=["sample", "amplicon"], inplace=True)
print(amp_cov)

scale = alt.Scale(type="linear", domain=[0,10,11,19,20,5000], range=["red", "red", "yellow", "yellow", "green", "green"])

plot = alt.Chart(amp_cov).mark_rect().encode(
    alt.X("amplicon:N", axis=alt.Axis(title="Amplicon number (n = 154)")),
    alt.Y("sample:O", axis=alt.Axis(title="Sampling location and week")),
    color=alt.Color('coverage:Q', scale=scale, legend=None)
).configure_axis(
    labelFontSize=20,
    titleFontSize=30
)

plot.save(snakemake.output.plot)