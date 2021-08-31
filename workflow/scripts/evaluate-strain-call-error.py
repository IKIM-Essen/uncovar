from os import name
import pandas as pd
import numpy as np


def extract_mixture_sample(path, prefix, separator, percentage, mix):
    path = path.split(prefix + separator)[-1]
    path = path.split(".strains")[0]
    path = path.replace("-", ".")
    path = path.split(separator)
    path = [ele.split(percentage) for ele in path]
    df = pd.DataFrame(path, columns=["target_id", "true_fraction"])
    df["mix"] = mix
    return df


def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


def eval_kallisto_error(
    kallisto_tsvs, sm_output, max_reads, prefix, separator, percentage
):
    results_df = pd.DataFrame()

    for i, path in enumerate(kallisto_tsvs):
        kallisto_df = pd.read_csv(path, delimiter="\t")
        kallisto_df["mix"] = i
        kallisto_df = kallisto_df.rename(columns={"fraction": "est_fraction"})

        org_mix_df = extract_mixture_sample(path, prefix, separator, percentage, i)
        org_mix_df["true_fraction"] = org_mix_df["true_fraction"].astype(int)
        org_mix_df["true_fraction"] = org_mix_df["true_fraction"] / 100
        org_mix_df["true_counts"] = round(org_mix_df["true_fraction"] * int(max_reads))

        kallisto_df = kallisto_df.merge(org_mix_df, how="outer").fillna(0)

        results_df = results_df.append(kallisto_df)

    for sample in results_df["mix"].unique():
        sample_rmse = rmse(
            results_df[results_df["mix"] == sample]["est_fraction"],
            results_df[results_df["mix"] == sample]["true_fraction"],
        )
        results_df.loc[results_df["mix"] == sample, "rmse"] = sample_rmse

    results_df.set_index(["mix", "rmse", "target_id"], inplace=True)
    results_df.to_csv(sm_output, sep="\t")


def eval_pangolin_error(
    pangolin_csvs, sm_output, max_reads, prefix, separator, percentage
):
    results_df = pd.DataFrame()

    for i, path in enumerate(pangolin_csvs):
        pangolin_df = pd.read_csv(path, delimiter=",")
        pangolin_df.rename(
            columns={"lineage": "target_id", "scorpio_support": "est_fraction"},
            inplace=True,
        )
        pangolin_df.drop(
            columns=["taxon", "pangoLEARN_version", "status", "note"], inplace=True
        )
        pangolin_df["mix"] = int(i)

        org_mix_df = extract_mixture_sample(path, prefix, separator, percentage, i)
        org_mix_df["true_fraction"] = org_mix_df["true_fraction"].astype(int)
        org_mix_df["true_fraction"] = org_mix_df["true_fraction"] / 100
        org_mix_df["true_counts"] = round(org_mix_df["true_fraction"] * int(max_reads))
        pangolin_df = pangolin_df.merge(org_mix_df, how="outer").fillna(0)

        results_df = results_df.append(pangolin_df, ignore_index=True)

    for sample in results_df["mix"].unique():
        sample_rmse = rmse(
            results_df[results_df["mix"] == sample]["est_fraction"],
            results_df[results_df["mix"] == sample]["true_fraction"],
        )
        results_df.loc[results_df["mix"] == sample, "rmse"] = sample_rmse

    results_df.set_index(["mix", "rmse", "target_id"], inplace=True)
    results_df.to_csv(sm_output, sep="\t")


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    max_reads = snakemake.params.get("max_reads", "")
    prefix = snakemake.params.get("prefix", "")
    separator = snakemake.params.get("separator", "")
    percentage = snakemake.params.get("percentage", "")

    if snakemake.wildcards.caller == "pangolin":
        eval_pangolin_error(
            snakemake.input,
            snakemake.output[0],
            max_reads,
            prefix,
            separator,
            percentage,
        )
    elif snakemake.wildcards.caller == "kallisto":
        eval_kallisto_error(
            snakemake.input,
            snakemake.output[0],
            max_reads,
            prefix,
            separator,
            percentage,
        )
