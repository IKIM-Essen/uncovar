sys.stderr = open(snakemake.log[0], "w")

from collections import Counter
from pysam import AlignmentFile


def get_softclip(record, idx):
    op, n = record.cigartuples[idx]
    if op == 4: # softclip
        return record.query_sequence[:idx * n]

softclipped_sequences = Counter()


# read sample BAM of reads mapping against Wuhan
with AlignmentFile(snakemake.input[0]) as bam:
    for record in bam:
        if record.is_unmapped:
            continue
        # extract softclipped sequences
        clipped_seq_start = get_softclip(record, 0)
        clipped_seq_end = get_softclip(record, -1)
        # we'll need to see if this is ok or too strict and whether we should count based on an inner substring or so
        softclipped_sequences[clipped_seq_start] += 1
        softclipped_sequences[clipped_seq_end] += 1

# print top 20 sequences
with open(snakemake.output[0], "w") as out_file:
    [print(f"{seq}: {count}", sep="\n", file=out_file) for seq, count in softclipped_sequences.most_common(20)]



# TODO additionally, it makes sense to also print the number and length of softclips in total
