import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

SUFFIXES = ["-medaka", "-nanopolish"]

happy_outputs = []

for happy_path, meta, platform in zip(
    snakemake.input, snakemake.params.metadata, snakemake.params.platforms
):
    workflow, sample = meta.split(",")

    happy = pd.read_csv(happy_path)

    happy["Mode"] = ""

    for suffix in SUFFIXES:
        if workflow.endswith(suffix):
            workflow = workflow.removesuffix(suffix)
            happy["Mode"] = suffix.removeprefix("-")

    happy["Workflow"] = workflow
    happy["Sample"] = sample
    happy["Platform"] = platform
    happy["Coverage"] = snakemake.wildcards.cov

    happy_outputs.append(happy)

happy_outputs = pd.concat(happy_outputs, ignore_index=True)
happy_outputs = happy_outputs[happy_outputs["Filter"] != "ALL"]

happy_outputs.to_csv(snakemake.output[0], sep="\t")
