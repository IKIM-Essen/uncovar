import subprocess as sp
import os
from pathlib import Path

import numpy as np
import pysam

from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

with pysam.AlignmentFile(snakemake.input.bam, "rb") as infile:
    a, c, g, t = infile.count_coverage()
    median_coverage = np.median(a + c + g + t)

with open(snakemake.output.depth) as outfile:
    print(median_coverage, file=outfile)

# see https://bitbucket.org/jbaaijens/savage/src/master/
patch_num = median_coverage / 1000

os.makedirs(snakemake.output[0])
indir = "../" * len(Path(snakemake.output[0]).parts)

if snakemake.params.denovo:
    ref = f"--ref {indir}{snakemake.input.ref}"
else:
    ref = ""

shell(
    "cd {snakemake.output.assembly}; savage -t {snakemake.threads} {ref} "
    "--split {patch_num} --p1 {indir}{snakemake.input.fastq1} --p2 {indir}{snakemake.input.fastq2} "
    "2> {log}"
)