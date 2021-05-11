import sys

sys.stderr = open(snakemake.log[0], "w")

import json
from os.path import join, isfile

import numpy as np
import pandas as pd


def extrace_strains_from_provision(
    path_to_provision: str, path_to_strain_summary: str, path_to_strain_genomes: str
):
    # select strain genomes
    provision = pd.DataFrame()
    chunks = pd.read_json(path_to_provision, lines=True, chunksize=10000)
    for i, chunk in enumerate(chunks):
        print(f"Parsing chunk {i}", file=sys.stderr)
        provision = provision.append(select_oldest_strains(chunk), ignore_index=True)
    provision = select_oldest_strains(provision)

    # save strain genomes
    provision["covv_lineage_fasta"] = provision["covv_lineage"].values + ".fasta"
    np.vectorize(write_sequence)(
        provision["covv_lineage"].values,
        provision["covv_lineage_fasta"].values,
        provision["sequence"].values,
        path_to_strain_genomes,
    )

    # save strain genome summary
    provision["covv_lineage_fasta"] = np.vectorize(join)(
        path_to_strain_genomes, provision["covv_lineage_fasta"].values
    )
    provision["covv_lineage_fasta"].to_csv(
        path_to_strain_summary, sep="\n", header=False, index=False
    )
    print(provision, file=sys.stderr)


def select_oldest_strains(df: pd.DataFrame):
    # covv_lineage -> pangolin lineage
    # n_content -> share of nan in seq
    # covv_subm_date -> submission date of seq
    # covv_host -> host of the seq
    # is_complete -> seq ist complete

    preselect_filter = (
        (df["covv_host"] == "Human")
        & (df["is_complete"] == True)
        & (df["covv_lineage"] != "None")
    )
    cols_of_interesst = ["covv_lineage", "n_content", "covv_subm_date"]
    df = df.copy()
    df.dropna(subset=cols_of_interesst, inplace=True)
    df = df[preselect_filter]
    df.sort_values(by=cols_of_interesst, inplace=True)
    df.drop_duplicates(subset=["covv_lineage"], inplace=True)
    return df


def write_sequence(
    covv_lineage: str, covv_lineage_fasta: str, sequence: str, out_path: str
):
    genome_file = join(out_path, covv_lineage_fasta)
    if not isfile(genome_file):
        with open(genome_file, "w") as f:
            f.write(f">{covv_lineage}\n")
            f.write(f"{sequence}\n")


if __name__ == "__main__":
    extrace_strains_from_provision(
        path_to_provision=snakemake.input[0],
        path_to_strain_summary=snakemake.output[0],
        path_to_strain_genomes=snakemake.params.get("save_strains_to", ""),
    )
