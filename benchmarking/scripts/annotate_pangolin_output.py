import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

pangolin_call = pd.read_csv(snakemake.input[0])

pangolin_call["Workflow"] = snakemake.wildcards.workflow
pangolin_call["Sample"] = snakemake.wildcards.sample
pangolin_call["Platform"] = snakemake.wildcards.tech

pangolin_call.to_csv(snakemake.output[0])
