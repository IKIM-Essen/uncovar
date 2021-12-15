# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import pysam

# No samples make it past the qulaity filter
if snakemake.input.contigs == "resources/genomes/main.fasta":
    # write empty fasta file
    with open(snakemake.output.fasta, "w") as outfile:
        pass

    # Write empty csv file
    pd.DataFrame(
        columns=[
            "SENDING_LAB",
            "DATE_DRAW",
            "SEQ_TYPE",
            "SEQ_REASON",
            "SAMPLE_TYPE",
            "PUBLICATION_STATUS",
            "OWN_FASTA_ID",
        ]
    ).to_csv(snakemake.output.table, sep=";", index=False)

else:
    # Aggregating fasta files
    sequence_names = []
    for file in snakemake.input.contigs:
        with pysam.FastxFile(file) as infile, open(
            snakemake.output.fasta, "a"
        ) as outfile:
            for entry in infile:
                print(f">{entry.name}", file=outfile)
                print(entry.sequence, file=outfile)
                sequence_names.append(entry.name)

    # Creating csv-table
    csv_table = pd.DataFrame(
        {
            "SENDING_LAB": 10259,
            "DATE_DRAW": "",
            "SEQ_TYPE": snakemake.params.seq_type,
            "SEQ_REASON": "N",
            "SAMPLE_TYPE": "s001",
            "PUBLICATION_STATUS": "N",
            "OWN_FASTA_ID": sequence_names,
        }
    )

    csv_table.sort_values(by="OWN_FASTA_ID", inplace=True)
    csv_table.to_csv(snakemake.output.table, sep=";", index=False)
