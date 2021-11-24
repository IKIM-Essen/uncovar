import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from pandas._typing import FilePathOrBuffer

summary = pd.DataFrame()


def register_quality_data(path_to_type_summary: FilePathOrBuffer, assembly_type: str):
    if path_to_type_summary != "resources/genomes/main.fasta":
        global summary
        quality_data = pd.read_csv(path_to_type_summary, sep="\t", index_col="Sample")
        quality_data["filter"] = (quality_data["identity"] > snakemake.params.min_identity) and (quality_data["n_share"] < snakemake.params.max_n)
        quality_data.rename(
            columns={
                "identity": "{}: Identity".format(assembly_type),
                "n_share": "{}: Share N".format(assembly_type),
                "filter": "{}: Pass Filter".format(assembly_type),
            },
            inplace=True,
        )
        summary = pd.concat([summary, quality_data], axis=1)


register_quality_data(snakemake.input.de_novo, "De Novo")
register_quality_data(snakemake.input.pseudo, "Pseudo")
register_quality_data(snakemake.input.consensus, "Consensus")

summary = summary.applymap(lambda x: "{:,.2f}%".format(x * 100))

summary.to_csv(snakemake.output[0])
