# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys

sys.stderr = open(snakemake.log[0], "w")
# parameter = snakemake.params.get("parameter", "")

from shutil import copyfile

from Bio import SeqIO


def is_fasta(filename):
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file


def check_contigs(sm_input, sm_output):
    if is_fasta(sm_input):
        copyfile(sm_input, sm_output)
    else:
        with open(sm_output, "w") as write_handle:
            write_handle.write(f">filler-contig\nN")


if __name__ == "__main__":
    check_contigs(snakemake.input[0], snakemake.output[0])
