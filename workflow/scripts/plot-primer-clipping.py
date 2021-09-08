import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import pysam
from intervaltree import IntervalTree
import altair as alt

# read primer bedpe to df
PRIMER = pd.read_csv(snakemake.params.get("bed", ""), delimiter="\t", header=None)
PRIMER.drop(PRIMER.columns[[0, 3]], axis=1, inplace=True)
PRIMER.columns = ["p1_start", "p1_end", "p2_start", "p2_end"]

# convert df to interval trees
primer_intervals = IntervalTree()
no_primer_intervals = IntervalTree()
for index, row in PRIMER.iterrows():
    primer_intervals[row["p1_start"] : row["p2_end"] + 1] = (
        row["p1_start"],
        row["p2_end"] + 1,
    )
    no_primer_intervals[row["p1_end"] + 1 : row["p2_start"]] = (
        row["p1_end"] + 1,
        row["p2_start"],
    )


def iter_with_samples(inputfiles):
    return zip(snakemake.params.samples, inputfiles)


def count_intervals(file):
    with pysam.AlignmentFile(file, "rb") as bam:
        counter_primer = 0
        counter_no_primer = 0
        counter_primer_within = 0
        counter_no_primer_within = 0
        counter_nothing = 0
        mate_pair_intervals = {}
        for read in bam.fetch():
            if not mate_pair_intervals.get(read.query_name):
                mate_pair_intervals[read.query_name] = [read.reference_start]
            else:
                mate_pair_intervals[read.query_name].append(read.reference_end)
        for pair in mate_pair_intervals:
            if (
                len(mate_pair_intervals[pair]) > 1
                and mate_pair_intervals[pair][0] != None
                and mate_pair_intervals[pair][1] != None
            ):
                if primer_intervals.envelop(
                    mate_pair_intervals[pair][0], mate_pair_intervals[pair][1] + 1
                ):
                    if (
                        sorted(
                            primer_intervals.envelop(
                                mate_pair_intervals[pair][0],
                                mate_pair_intervals[pair][1] + 1,
                            )
                        )[0].begin
                        == mate_pair_intervals[pair][0]
                        and sorted(
                            primer_intervals.envelop(
                                mate_pair_intervals[pair][0],
                                mate_pair_intervals[pair][1] + 1,
                            )
                        )[0].end
                        == mate_pair_intervals[pair][1] + 1
                    ):
                        counter_primer += 1
                    else:
                        counter_primer_within += 1
                elif no_primer_intervals.envelop(
                    mate_pair_intervals[pair][0] + 1, mate_pair_intervals[pair][1]
                ):
                    if (
                        sorted(
                            no_primer_intervals.envelop(
                                mate_pair_intervals[pair][0] + 1,
                                mate_pair_intervals[pair][1],
                            )
                        )[0].begin
                        == mate_pair_intervals[pair][0] + 1
                        and sorted(
                            no_primer_intervals.envelop(
                                mate_pair_intervals[pair][0] + 1,
                                mate_pair_intervals[pair][1],
                            )
                        )[0].end
                        == mate_pair_intervals[pair][1]
                    ):
                        counter_no_primer += 1
                    else:
                        counter_no_primer_within += 1
                else:
                    counter_nothing += 1
            else:
                counter_nothing += 1
        counters = pd.DataFrame(
            {
                "n_count": [
                    counter_primer,
                    counter_primer_within,
                    counter_no_primer,
                    counter_no_primer_within,
                    counter_nothing,
                ],
                "class": [
                    "uncut primer exact",
                    "uncut primer within",
                    "cut primer exact",
                    "cut primer within",
                    "no mathing win",
                ],
            }
        )
        return counters


def plot_classes(counters):
    bars = (
        alt.Chart(counters)
        .mark_bar()
        .encode(
            y="class",
            x="n_count",
            row=alt.Row("sample:N"),
            column=alt.Column("state:N", sort="descending"),
        )
    )
    text = bars.mark_text(
        align="left",
        baseline="middle",
        dx=3,  # Nudges text to right so it doesn't appear on top of the bar
    ).encode(text="n_count", row=alt.Row("sample:N"), column=alt.Column("state:N"))
    return bars, text


all_df = pd.DataFrame()
for sample, file in iter_with_samples(snakemake.input.unclipped):
    counts_before = count_intervals(file)
    counts_before["sample"] = sample
    counts_before["state"] = "before"
    all_df = all_df.append(counts_before, ignore_index=True)

for sample, file in iter_with_samples(snakemake.input.clipped):
    counts_after = count_intervals(file)
    counts_after["sample"] = sample
    counts_after["state"] = "after"
    all_df = all_df.append(counts_after, ignore_index=True)

bars, text = plot_classes(all_df)

(bars).properties(title="Amplicon matching").save(snakemake.output.plot)
