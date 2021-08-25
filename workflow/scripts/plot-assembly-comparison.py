import sys
from altair.vegalite.v4.schema.channels import Opacity
import pandas as pd
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

def register_contig_lengths(assemblies_init, assemblies_final, N50, assembler, data):
    for sample, amplicon_state, file1, file2, quast in zip(snakemake.params.samples, snakemake.params.amplicon_state, assemblies_init, assemblies_final, N50):
        with pysam.FastxFile(file1) as infile_init:
            with pysam.FastxFile(file2) as infile_final:
                quastDf =  pd.read_csv(quast,  sep='\t')
                data = data.append({
                    'Sample': sample,
                    'Assembler': assembler,
                    'Amplicon':amplicon_state, 
                    'bp': max(len(contig.sequence) for contig in infile_init),
                    'Contig': 'initial',
                    'N50': quastDf.loc[0, 'N50'],
                    'Genome fraction (%)': quastDf.loc[0, 'Genome fraction (%)'] if 'Genome fraction (%)' in quastDf.columns else float('nan'),
                    },
                    ignore_index=True
                )
                data = data.append({
                    'Sample': sample,
                    'Assembler': assembler,
                    'Amplicon':amplicon_state, 
                    'bp': max(len(contig.sequence) for contig in infile_final),
                    'Contig': 'final',
                    },
                    ignore_index=True
                )
    return data

dict_init = make_assembler_lists(snakemake.input.initial)
dict_final = make_assembler_lists(snakemake.input.final)
dict_N50 = make_assembler_lists(snakemake.input.quast)
for assembler in dict_init:
    data=register_contig_lengths(dict_init[assembler], dict_final[assembler], dict_N50[assembler], assembler, data)

data.replace({'Amplicon': {0: 'Shotgun', 1: 'Amplicon'}}, inplace=True)
data.replace({'Assembler': {
    "megahit-std": 'MEGAHIT',
    "megahit-meta-large": "MEGAHIT meta-large",
    "megahit-meta-sensitive": "MEGAHIT meta-sensitive",
    "trinity": "Trinty",
    "velvet": "Velvet",
    "metaspades": "SPAdes meta",
    "spades": "SPAdes",
    "coronaspades": "SPAdes corona",
    "rnaviralspades": "SPAdes RNA viral",
     }}, inplace=True
)

data.to_csv(snakemake.output[1])

height=300
width=50

plot_bp = alt.Chart(data).encode(
    y=alt.Y('bp',  scale=alt.Scale(domain=[0, 35000], clamp=True), axis=alt.Axis(tickCount=9)),
    x=alt.X('Contig', title=None, sort='descending'),
    color=alt.Color('Assembler', scale=alt.Scale(scheme='turbo'), legend=None),
).properties(height=height, width=width)

combined_bp = plot_bp.mark_boxplot(opacity=0.5) + plot_bp.mark_point(opacity=0.5, filled=True)
combined_bp = combined_bp.facet(
    column= alt.Column('Assembler:N',
        title="",
        header=alt.Header(labelAngle=-45, labelOrient='bottom', labelPadding=-5)
    ), 
    row=alt.Row('Amplicon:N',
        title="",
        sort='descending',
    ),
).configure_axis(
    labelFontSize=12,
    titleFontSize=12
).configure_view(
    stroke=None
)

plot_N50 = alt.Chart(data).mark_boxplot(opacity=0.5).encode(
    y=alt.Y('N50',  scale=alt.Scale(domain=[0, 35000], clamp=True), axis=alt.Axis(tickCount=9)),
    color=alt.Color('Assembler', scale=alt.Scale(scheme='turbo'), legend=None),
    column=alt.Column('Assembler:N',
        title="",
        header=alt.Header(labelAngle=-45, labelOrient='bottom', labelPadding=-5)
    ),
    row=alt.Row('Amplicon:N',
        title="",
        sort='descending',
    ), 
).configure_axis(
    grid=False,
    labelFontSize=12,
    titleFontSize=12
).configure_view(
    stroke=None
).properties(
    height=height, width=width
)

plot_genome_frac = alt.Chart(data).mark_point(opacity=0.5).encode(
    x=alt.X('jitter:Q',
        title=None,
        axis=alt.Axis(ticks=True, grid=False, labels=False),
        scale=alt.Scale(),
    ),
    y=alt.Y('Genome fraction (%)',  scale=alt.Scale(domain=[0, 100]), axis=alt.Axis(tickCount=9)),
    color=alt.Color('Assembler', scale=alt.Scale(scheme='turbo'), legend=None),
    column=alt.Column('Assembler:N',
        title="",
        header=alt.Header(labelAngle=-45, labelOrient='bottom', labelPadding=-5)
    ),
    row=alt.Row('Amplicon:N',
        title="",
        sort='descending',
    ), 
).configure_axis(
    grid=False,
    labelFontSize=12,
    titleFontSize=12
).transform_calculate(
    jitter='sqrt(-2*log(random()))*cos(2*PI*random())'
).configure_view(
    stroke=None
).configure_axis(
    labelFontSize=12,
    titleFontSize=12
).properties(height=height, width=width)

combined_bp.save(snakemake.output[0])
plot_N50.save(snakemake.output[2])
plot_genome_frac.save(snakemake.output[3])