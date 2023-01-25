import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

calls = []
for path in snakemake.input:
    calls.append(pd.read_csv(path, sep="\t"))

calls = pd.concat(calls, ignore_index=True)
calls.to_csv(snakemake.output[0], sep="\t", index=False)
