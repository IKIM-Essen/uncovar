import sys

sys.stderr = open(snakemake.log[0], "w")

from os import name

import numpy as np
import pandas as pd


def extract_mixture_sample(path, prefix, separator, percentage, mix, max_reads):
    path = path.split(prefix + separator)[-1]
    path = path.split(".strains")[0]
    path = path.replace("-", ".")
    path = path.split(separator)
    path = [ele.split(percentage) for ele in path]
    df = pd.DataFrame(path, columns=["target_id", "true_fraction"])
    df["mix"] = mix
    df["true_fraction"] = df["true_fraction"].astype(int)
    df["true_fraction"] = df["true_fraction"] / 100
    df["true_counts"] = round(df["true_fraction"] * int(max_reads))
    return df


def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


def load_kallisto_df(i, path):
    kallisto_df = pd.read_csv(path, delimiter="\t")
    kallisto_df["mix"] = i
    kallisto_df = kallisto_df.rename(columns={"fraction": "est_fraction"})
    return kallisto_df


def load_pangolin_df(i, path):
    pangolin_df = pd.read_csv(path, delimiter=",")
    pangolin_df.rename(
        columns={"lineage": "target_id", "scorpio_support": "est_fraction"},
        inplace=True,
    )
    pangolin_df.drop(
        columns=["taxon", "pangoLEARN_version", "status", "note"], inplace=True
    )
    pangolin_df["mix"] = i
    return pangolin_df


def eval_error(paths, sm_output, max_reads, prefix, separator, percentage, load_func):
    results_df = pd.DataFrame()

    for i, path in enumerate(paths):
        df = load_func(i, path)

        org_mix_df = extract_mixture_sample(
            path, prefix, separator, percentage, i, max_reads
        )

        df = df.merge(org_mix_df, how="outer").fillna(0)

        results_df = results_df.append(df)

    for sample in results_df["mix"].unique():
        sample_rmse = rmse(
            results_df[results_df["mix"] == sample]["est_fraction"],
            results_df[results_df["mix"] == sample]["true_fraction"],
        )
        results_df.loc[results_df["mix"] == sample, "rmse"] = sample_rmse

    results_df.set_index(["mix", "rmse", "target_id"], inplace=True)
    results_df.to_csv(sm_output, sep="\t")


max_reads = snakemake.params.max_reads
prefix = snakemake.params.prefix
separator = snakemake.params.separator
percentage = snakemake.params.percentage

if snakemake.wildcards.caller == "pangolin":
    load_func = load_pangolin_df
elif snakemake.wildcards.caller == "kallisto":
    load_func = load_kallisto_df
else:
    raise ValueError("unexpected value for wildcards.caller")

eval_error(
    snakemake.input,
    snakemake.output[0],
    max_reads,
    prefix,
    separator,
    percentage,
    load_func,
)
