import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import pysam

number_of_bases = 0

for path in snakemake.input:
    with pysam.FastxFile(path) as fh:
        for entry in fh:
            number_of_bases += len(entry.sequence)

number = pd.DataFrame(
    {"Sample": [snakemake.wildcards.sample], "# Bases": [number_of_bases]}
).to_csv(snakemake.output[0], sep="\t", index=False)
