import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd


pangolin_calls = []

assert len(snakemake.input) == len(snakemake.params.samples)
assert len(snakemake.input) == len(snakemake.params.stages)

for path, sample, stage in zip(
    snakemake.input, snakemake.params.samples, snakemake.params.stages
):
    pangolin_call = pd.read_csv(path)
    pangolin_call["sample"] = sample
    pangolin_call["stage"] = stage

    pangolin_calls.append(pangolin_call)

pangolin_calls_by_stage = pd.concat(pangolin_calls, axis=0, ignore_index=True)
pangolin_calls_by_stage = pangolin_calls_by_stage.pivot(
    index="sample", columns="stage", values="lineage"
).reset_index()

only_pseudo_assembly = ["scaffold", "polished", "masked-polished", "pseudo"]
only_consensus_assembly = [
    "scaffold",
    "polished",
    "masked-polished",
    "consensus",
    "masked-consensus",
]
both_fallbacks = list(
    set(only_pseudo_assembly).symmetric_difference(set(only_consensus_assembly))
)

columns = pangolin_calls_by_stage.columns.to_list()

if all(stage in columns for stage in both_fallbacks):
    sorted_columns = [
        "sample",
        "scaffold",
        "polished",
        "masked-polished",
        "consensus",
        "masked-consensus",
        "pseudo",
    ]
elif all(stage in columns for stage in only_pseudo_assembly):
    sorted_columns = ["sample", "scaffold", "polished", "masked-polished", "pseudo"]
elif all(stage in columns for stage in only_consensus_assembly):
    sorted_columns = [
        "sample",
        "scaffold",
        "polished",
        "masked-polished",
        "consensus",
        "masked-consensus",
    ]

pangolin_calls_by_stage = pangolin_calls_by_stage[sorted_columns]

pangolin_calls_by_stage.rename(
    columns={
        "sample": "Sample",
        "scaffold": "Scaffolded Seq",
        "polished": "Polished Seq",
        "masked-polished": "Masked Polished Seq",
        "pseudo": "Pseudo Seq",
        "consensus": "Consensus Seq",
        "masked-consensus": "Masked Consensus Seq",
    },
    inplace=True,
)

pangolin_calls_by_stage.fillna("-", inplace=True)
pangolin_calls_by_stage.replace({"None": "Call failed"}, inplace=True)
pangolin_calls_by_stage.sort_values(by="Sample", inplace=True)

print(pangolin_calls_by_stage)

pangolin_calls_by_stage.to_csv(snakemake.output[0], index=False)
