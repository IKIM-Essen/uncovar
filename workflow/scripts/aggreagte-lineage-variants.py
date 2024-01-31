# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys
import pandas as pd
import numpy as np
import gffutils
from natsort import natsort_keygen

sys.stderr = open(snakemake.log[0], "w")

data = pd.DataFrame()

for file, sample in zip(snakemake.input.csv, snakemake.params.sample):
    data = pd.concat(
        [
            data,
            pd.read_csv(file, usecols=["Mutations","Probability","Frequency","ReadDepth"])
        ]
    )
    data.rename(columns={"Probability": sample + ": " + "Probability"}, inplace=True)
    data.rename(columns={"Frequency": sample + ": " + "VAF"}, inplace=True)
    data.rename(columns={"ReadDepth": sample + ": " + "Read Depth"}, inplace=True)

data = data.groupby(["Mutations"]).agg(func={column: np.max for column in data.columns})
# add feature column for sorting
data["Features"] = data.index.to_series().str.extract(r"(.+)[:].+|\*")

# position of variant for sorting and change type
data["Position"] = data.index.to_series().str.extract(
    r"(.*:?[A-Z]+|\*$|-)([0-9]+)([A-Z]+$|\*$|-)$"
)[1]
# data = data.astype({"Position": "int64"})

# generate sorting list from .gff with correct order of features
gff = gffutils.create_db(snakemake.input.annotation, dbfn=":memory:")
gene_start = {gene["gene_name"][0]: gene.start for gene in gff.features_of_type("gene")}
sorter = [k[0] for k in sorted(gene_start.items(), key=lambda item: item[1])]
sorterIndex = dict(zip(sorter, range(len(sorter))))
data["Features_Rank"] = data["Features"].map(sorterIndex)

data["Features_Rank"].replace(np.NaN, 12, inplace=True)

data.drop(index=["Lineage", "Similarity"], inplace=True)
data.sort_values(
    by=["Features_Rank", "Position"],
    ascending=[True, True],
    na_position="last",
    inplace=True,
    key=natsort_keygen(),
)

data.drop(columns=["Mutations", "Features_Rank"], index=[], inplace=True)
print(data)
data.to_csv(snakemake.output[0])