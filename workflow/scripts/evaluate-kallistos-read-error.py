import pandas as pd

def extract_mixture_sample(path):
    path = path.split('mixture-sample-#')[-1]
    path = path.split('.strains')[0]
    path = path.replace("-", ".")
    path = path.split("#")
    path = [ele.split("=") for ele in path]
    return pd.DataFrame(path, columns=["target_id", "true_fraction"])


def eval_error(kallisto_tsvs, sm_output, max_reads):
    results_df = pd.DataFrame()
    
    for i, path in enumerate(kallisto_tsvs):
        kallisto_df = pd.read_csv(path, delimiter="\t")
        kallisto_df["sample"] = i
        kallisto_df.rename(columns={"fraction" : "est_fraction"})

        org_mix_df = extract_mixture_sample(path)
        org_mix_df["true_fraction"] =  org_mix_df["true_fraction"].astype(int)
        org_mix_df["true_fraction"] = org_mix_df["true_fraction"] / 100
        org_mix_df["true_counts"] = round(org_mix_df["true_fraction"]  * int(max_reads))

        kallisto_df = kallisto_df.merge(org_mix_df, how='outer').fillna(0)

        results_df=results_df.append(kallisto_df)
    
    results_df.set_index(["sample", "target_id"], inplace=True)
    print(results_df)


if __name__ == "__main__":
    sys.stderr = open(snakemake.log[0], "w")
    max_reads = snakemake.params.get("max_reads", "")

    eval_error(snakemake.input, snakemake.output[0], max_reads)
