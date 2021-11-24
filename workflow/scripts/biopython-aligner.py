import os
from numpy import diff
import pandas as pd
from Bio import SeqIO
from Bio import Align


def alignment(file1, file2):
    data = pd.DataFrame()
    seq1 = SeqIO.read(file1, "fasta")
    seq2 = SeqIO.read(file2, "fasta")
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignments = aligner.align(seq1.seq, seq2.seq)
    for alignment in alignments:
        print(alignment.shape[0])
        with open(snakemake.output[1], "w") as outfile:
            print(format(alignment, "sam"), file=outfile)
        aln = format(alignment, "psl").split("\t")
        data = data.append(
            {
                "Sample": str(snakemake.wildcards.sample),
                "Region": snakemake.wildcards.region,
                "Length": alignment.shape[1],
                "Matches(#)": int(aln[0]),
                "Matches(%)": int(aln[0]) / alignment.shape[1] * 100,
                "Matches w/o N's(%)": int(aln[0])
                / (alignment.shape[1] - int(aln[3]))
                * 100,
                "Mismatches(#)": int(aln[1]),
                "Mismatches(%)": int(aln[1]) / alignment.shape[1] * 100,
                "N(#)": int(aln[3]),
                "N(%)": int(aln[1]) / alignment.shape[1] * 100,
                "Gaps(#)": int(aln[4]),
                "Gaps(%)": int(aln[4]) / alignment.shape[1] * 100,
                "Aln Start(t)": int(aln[15]),
                "Aln End(t)": int(aln[16]),
            },
            ignore_index=True,
        )
        # print(data)
        return data


sanger_vs_genome = alignment(snakemake.input[0], snakemake.input[1])
sanger_vs_genome = sanger_vs_genome.set_index("Sample")
sanger_vs_genome.to_csv(snakemake.output[0])
