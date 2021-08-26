
import sys
from typing import Sequence
sys.stderr = open(snakemake.log[0], "w")

from Bio import SeqIO

def get_largest_contig(sm_input, sm_output):
    contigs = []

    for seq_record in SeqIO.parse(sm_input, "fasta"):
        contigs.append([seq_record.id, str(seq_record.seq), len(seq_record)])

    contigs.sort(key=lambda x: x[2])

    name = f">{contigs[-1][0]}"
    sequence = contigs[-1][1]

    with open(sm_output, "w") as writer:
        writer.write(
           f"{name}\n{sequence}"
        )

if __name__ == "__main__":
    get_largest_contig(snakemake.input[0], snakemake.output[0])