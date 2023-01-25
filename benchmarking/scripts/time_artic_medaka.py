import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

minion = pd.read_csv(snakemake.input.minion, sep="\t")
minion["source"] = "minion"

guppyplex = pd.read_csv(snakemake.input.guppyplex, sep="\t")
guppyplex["source"] = "guppyplex"

pd.concat([minion, guppyplex]).to_csv(snakemake.output[0], sep="\t", index=False)
