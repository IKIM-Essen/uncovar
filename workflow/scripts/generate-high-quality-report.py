# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import pysam

# No samples make it past the quality filter
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
    include_flag = []
    sample_dict = {}
    for sample in snakemake.params.includeflag:
        sample_dict.update(sample)

    with open(snakemake.output.fasta, "w") as outfile:
        for file in snakemake.input.contigs:
            with pysam.FastxFile(file) as infile:
                for entry in infile:
                    sequence_names.append(entry.name)
                    to_include = int(sample_dict.get(entry.name))
                    include_flag.append(to_include)
                    if to_include:
                        print(f">{entry.name}", file=outfile)
                        print(entry.sequence, file=outfile)

    # Creating csv-table
    csv_table = pd.DataFrame(
        {
            "SENDING_LAB": snakemake.params.sending_lab_number,
            "DATE_DRAW": snakemake.params.date_draw,
            "SEQ_TYPE": snakemake.params.seq_type,
            "SEQ_REASON": "N",
            "SAMPLE_TYPE": "s001",
            "PUBLICATION_STATUS": "N",
            "OWN_FASTA_ID": sequence_names,
            "include": include_flag,
        }
    )

    # Only include samples with include flag
    csv_table["include"] = csv_table["include"].astype(int)
    csv_table = csv_table[csv_table["include"] == 1]
    csv_table.drop(columns=["include"], inplace=True)

    # Final touches
    csv_table.sort_values(by="OWN_FASTA_ID", inplace=True)
    csv_table.to_csv(snakemake.output.table, sep=";", index=False)
