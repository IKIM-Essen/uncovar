# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

from snakemake.shell import shell
import pysam

exclude = (
    "-x {}".format(snakemake.input.exclude)
    if snakemake.input.get("exclude", "")
    else ""
)

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with pysam.AlignmentFile(snakemake.input.sample) as bam:
    if bam.mapped < 10000:
        # Not enough reads to perform SV calling.
        # Output empty BCF.
        header = pysam.VariantHeader()

        # Retrieve and record reference lengths.
        ref = pysam.FastaFile(snakemake.input.ref)
        for contig in ref.references:
            n = ref.get_reference_length(contig)
            header.add_line("##contig=<ID={},length={}>".format(contig, n))

        # Write BCF.
        with pysam.VariantFile(snakemake.output[0], "wb", header=header) as bcf:
            exit(0)


shell(
    "OMP_NUM_THREADS={snakemake.threads} delly call {extra} "
    "{exclude} -g {snakemake.input.ref} "
    "-o {snakemake.output[0]} {snakemake.input.sample} {log}"
)
