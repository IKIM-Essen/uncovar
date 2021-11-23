import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from pandas._typing import FilePathOrBuffer



summary = pd.DataFrame(index=snakemake.params.samples)

print(
    summary
)

def register_quality_data(path_to_type_summary: FilePathOrBuffer, assembly_type: str):
    if path_to_type_summary != "resources/genomes/main.fasta":
        quality_data = pd.read_csv(path_to_type_summary, sep="\t")
        quality_data.rename(columns={
            "identity" : "%s: Identity to Reference (%)".format(assembly_type),
            "n_share" : "%s: Share N in sequence (%)".format(assembly_type),
        }, inplace=True)
