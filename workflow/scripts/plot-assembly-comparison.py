import os
import sys
import pandas as pd
from pandas.core.indexes.base import Index
import pysam
import altair as alt

sys.stderr = open(snakemake.log[0], "w")

data = pd.DataFrame()

def make_assembler_lists(assemblies):
    dict = {}
    for assembler in snakemake.params.assembler:
        for file in assemblies:
            if assembler in file:
                if dict.get(assembler):
                    dict[assembler].append(file)
                else:
                    dict[assembler] = [file]
    return dict

def register_contig_lengths(assemblies_init, assemblies_final, assembler, data):
    for sample, file1, file2 in zip(snakemake.params.samples, assemblies_init, assemblies_final):
        with pysam.FastxFile(file1) as infile_init:
            with pysam.FastxFile(file2) as infile_final:
                data = data.append({
                    'Sample': sample,
                    'Assembler': assembler, 
                    'bp': max(len(contig.sequence) for contig in infile_init),
                    'Contig': 'initial',
                    },
                    ignore_index=True
                )
                data = data.append({
                    'Sample': sample,
                    'Assembler': assembler, 
                    'bp': max(len(contig.sequence) for contig in infile_final),
                    'Contig': 'final',
                    },
                    ignore_index=True
                )
    return data
                # data.loc[sample, "Assembler"] = assembler
                # data.loc[sample, "Initial Contig (bp)"] = max(len(contig.sequence) for contig in infile_init)
                # data.loc[sample, "Final Contig (bp)"] = max(len(contig.sequence) for contig in infile_final)

dict_init = make_assembler_lists(snakemake.input.initial)
dict_final = make_assembler_lists(snakemake.input.final)
for assembler in dict_init:
    data=register_contig_lengths(dict_init[assembler], dict_final[assembler], assembler, data)

# data.reset_index(inplace=True)


# init = alt.Chart(data).mark_boxplot(color='#57A44C').encode(
#     alt.Y('Initial Contig (bp)',
#           axis=alt.Axis(title='# bp')),
# )
plot1 = alt.Chart(data).mark_boxplot(color='#08A4C2').encode(
    alt.Y('bp'),
    color='Contig',
    x=('Assembler')
)

plot2 = alt.Chart(data).mark_boxplot().encode(
    x=alt.X('bp',  scale=alt.Scale(domain=[0, 30000]), axis=alt.Axis(tickCount=9)),
    y=alt.Y('Contig', title=None, sort='descending'),
    color=alt.Color('Contig', scale=alt.Scale(scheme='tableau20'), legend=None),
    row='Assembler',
    
)

plot3 =  alt.Chart(data).mark_circle(size=8).encode(
    x=alt.X(
        'bp',
        title=None,
        axis=alt.Axis(values=[0], ticks=True, grid=False, labels=False),
        scale=alt.Scale(),
    ),
    y=alt.Y('Contig'),
    color=alt.Color('Contig', legend=None),
    row=alt.Row(
        'Assembler',
        # header=alt.Header(
        #     labelAngle=-90,
        #     titleOrient='top',
        #     labelOrient='bottom',
        #     labelAlign='right',
        #     labelPadding=3,
        # ),
    ),
).transform_calculate(
    # Generate Gaussian jitter with a Box-Muller transform
    jitter='sqrt(-2*log(random()))*cos(2*PI*random())'
# ).configure_facet(
#     spacing=0
).configure_view(
    stroke=None
)

print(data)
# line = base.mark_line(stroke='#5276A7', interpolate='monotone').encode(
#     alt.Y('Initial Contig (bp) megahit',
#           axis=alt.Axis(title='Precipitation (inches)', titleColor='#5276A7'))

plot1.save(snakemake.output[0])
plot2.save(snakemake.output[0])
# (plot2+plot3).save(snakemake.output[0])
# with open(snakemake.output[0], "w") as outfile:

#     pass
