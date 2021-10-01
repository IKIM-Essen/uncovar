import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import pysam

# Aggregating fasta files
sequence_names = []
for file in snakemake.input.contigs:
    with pysam.FastxFile(file) as infile, open(snakemake.output.fasta, "a") as outfile:
        for entry in infile:
            print(entry.name, file=outfile)
            print(entry.sequence, file=outfile)
            sequence_names.append(entry.name)

# Creating csv-table
csv_table = pd.DataFrame(
    {
        "SENDING_LAB": 10259,
        "DATE_DRAW": "",
        "SEQ_TYPE": "ILLUMINA",
        "SEQ_REASON": "N",
        "SAMPLE_TYPE": "s001",
        "PUBLICATION_STATUS": "N",
        "OWN_FASTA_ID": sequence_names,
    }
)
csv_table.to_csv(snakemake.output.table, sep=";", index=False)
