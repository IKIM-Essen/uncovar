import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd


def aggregate_calls(sm_input, sm_output):
    all_calls = []
    for i, file in enumerate(sm_input):
        call = pd.read_csv(file, sep="\t")
        call["sample"] = i
        all_calls.append(call)

    all_calls = pd.concat(all_calls, axis=0, ignore_index=True)
    all_calls.to_csv(sm_output, sep="\t", index=False)


if __name__ == "__main__":
    aggregate_calls(snakemake.input, snakemake.output[0])
