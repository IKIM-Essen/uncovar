import sys

sys.stderr = open(snakemake.log[0], "w")

import re

import pandas as pd


def set_id(info_string: str, format: str) -> str:
    if "ID=" in info_string:
        return info_string

    infos = dict(x.split("=") for x in info_string.split(";"))

    info_string = f"ID={format}-{infos['Name']};" + info_string
    return info_string


def fix_transcipt(info_string: str, type: str, to_change: str, true_parent: str) -> str:
    if type != to_change:
        return info_string

    infos = dict(x.split("=") for x in info_string.split(";"))

    return re.sub("Parent=?.*;", f'Parent={true_parent}-{infos["Name"]};', info_string)


gff = pd.read_csv(snakemake.input[0], sep="\t", header=None)

info_column = gff.columns[-1]
format_column = gff.columns[2]

exons = gff.loc[gff[format_column] == "CDS"].copy()
exons.replace(to_replace="CDS", value="exon", inplace=True)
gff = pd.concat([gff, exons])


gff[info_column] = gff.apply(lambda x: set_id(x[info_column], x[format_column]), axis=1)
gff[info_column] = gff.apply(
    lambda x: fix_transcipt(
        x[info_column], x[format_column], to_change="CDS", true_parent="transcript"
    ),
    axis=1,
)
gff[info_column] = gff.apply(
    lambda x: fix_transcipt(
        x[info_column], x[format_column], to_change="exon", true_parent="transcript"
    ),
    axis=1,
)


gff.loc[gff[format_column] == "gene", info_column] = gff.loc[
    gff[format_column] == "gene", info_column
].apply(lambda x: x + ";biotype=protein_coding")

gff[format_column].replace(to_replace="transcript", value="mRNA", inplace=True)

gff.to_csv(snakemake.output[0], header=False, index=False, sep="\t")
