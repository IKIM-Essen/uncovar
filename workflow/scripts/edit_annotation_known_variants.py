import pandas as pd

input_gff = pd.read_csv(
    snakemake.input[0],
    delimiter="\t",
    names=[
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ],
)

input_gff["seqname"] = "MN908947.3"
input_gff["source"] = "nextclade"

CDS_gff = input_gff.copy()
CDS_gff["feature"] = "CDS"
CDS_gff["frame"] = "0"

exon_gff = input_gff.copy()
exon_gff["feature"] = "exon"

mRNA_gff = input_gff.copy()
mRNA_gff["feature"] = "mRNA"
mRNA_gff["attribute"] += ";biotype=protein_coding"

output_df = pd.concat([input_gff, CDS_gff, exon_gff, mRNA_gff])
output_df.sort_values(by=["start", "feature"], inplace=True)
print(output_df)

output_df.to_csv(snakemake.output[0], sep="\t", index=False, header=False)
