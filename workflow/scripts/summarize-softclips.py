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
        softclipped_sequences[clipped_seq_start] += 1
        softclipped_sequences[clipped_seq_end] += 1

# print top 20 sequences
print(*map("{}: {}".format, softclipped_sequences.most_common(20)), sep="\n", file=snakemake.output[0])
