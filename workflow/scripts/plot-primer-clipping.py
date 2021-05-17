from numpy.core.numeric import base_repr
import pandas as pd
import pysam
from intervaltree import Interval, IntervalTree
import altair as alt

# read primer bedpe to df
PRIMER = pd.read_csv("/local/data/alex/testing/snakemake-workflow-sars-cov2/resources/primer.bedpe", delimiter="\t", header=None)
PRIMER.drop(PRIMER.columns[[0, 3]], axis=1, inplace=True)
PRIMER.columns = ['p1_start', 'p1_end', 'p2_start', 'p2_end']

# convert df to interval trees
primer_intervals = IntervalTree()
no_primer_intervals = IntervalTree()
for index, row in PRIMER.iterrows():
    primer_intervals[row['p1_start']: row['p2_end']+1] = (row['p1_start'], row['p2_end']+1)
    no_primer_intervals[row['p1_end']+1: row['p2_start']] = (row['p1_end']+1, row['p2_start'])


def count_intervals(file, state):
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
            if mate_pair_intervals[pair][0] != None and mate_pair_intervals[pair][1] != None:
                if primer_intervals.envelop(mate_pair_intervals[pair][0], mate_pair_intervals[pair][1]+1):
                    if sorted(primer_intervals.envelop(mate_pair_intervals[pair][0], mate_pair_intervals[pair][1]+1))[0].begin == mate_pair_intervals[pair][0] \
                    and sorted(primer_intervals.envelop(mate_pair_intervals[pair][0], mate_pair_intervals[pair][1]+1))[0].end == mate_pair_intervals[pair][1]+1:
                        counter_primer += 1
                    else:
                        counter_primer_within += 1
                elif no_primer_intervals.envelop(mate_pair_intervals[pair][0]+1, mate_pair_intervals[pair][1]):
                    if sorted(no_primer_intervals.envelop(mate_pair_intervals[pair][0]+1, mate_pair_intervals[pair][1]))[0].begin == mate_pair_intervals[pair][0]+1 \
                    and sorted(no_primer_intervals.envelop(mate_pair_intervals[pair][0]+1, mate_pair_intervals[pair][1]))[0].end == mate_pair_intervals[pair][1]:
                        counter_no_primer += 1
                    else:
                        counter_no_primer_within += 1
                else:
                    counter_nothing += 1
            else:
                counter_nothing += 1
        counters = pd.DataFrame({'n_count': [counter_primer, counter_primer_within, counter_no_primer, counter_no_primer_within, counter_nothing], state: ["w/ primer", "w/ primer within", "w/o primer", "w/o primer within", "none"]})
        return counters


def plot_classes(counters, state):
    bars = alt.Chart(counters).mark_bar().encode(
    # y=alt.Y("class", title=""),
    y=state,
    x='n_count',)
    text = bars.mark_text(
    align='left',
    baseline='middle',
    dx=3  # Nudges text to right so it doesn't appear on top of the bar
    ).encode(
    text='n_count'
    ) 
    return bars, text


counts_before = count_intervals("/local/data/alex/testing/snakemake-workflow-sars-cov2/results/2021-04-30/clipped-reads/57815.bam", "before")
bars_before, text_before = plot_classes(counts_before, "before")
counts_after = count_intervals("/local/data/alex/testing/snakemake-workflow-sars-cov2/results/2021-04-30/clipped-reads/57815.primerclipped.hard.bam", "after")
bars_after, text_after = plot_classes(counts_after, "after")

(bars_before + text_before & bars_after + text_after).properties(title='Amplicon matching').save("test.svg")