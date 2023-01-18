# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

# import required packages
import pandas as pd
import numpy as np
import re
import textwrap
import sys

sys.stderr = open(snakemake.log[0], "w")

# list of mutations effectiv in escaping the antibodies.
escaping_mutations = pd.read_json(snakemake.input.escaping_mutations).set_index("mAbs")
all_escaping_mutations = [
    item for sublist in escaping_mutations["Mutation"].tolist() for item in sublist
]

# list of factors
factors = pd.read_csv(snakemake.input.factors)
factors = factors[["Mutation", "S309", "COV2-2130", "COV2-2196"]].set_index("Mutation")

# list of mutations in samples
df = pd.read_csv(snakemake.input.allmutationsfound)
# filter for tests
# TODO keep all filenames ?list? # perhabs do list of sample ids and append missing ones
df = df[
    (df["Gene"] == "S")
    & (df["Frequency"] > 0)
    & (df["ReadDepth"] > 10)
    & ((df["ReadDepth"] * df["Frequency"]) > 1)
]

df = df.drop(columns=["Position", "ReadDepth", "Probability"])
df = df[["Sample_id", "Signature", "Gene", "Frequency", "Lineage"]].reset_index(
    drop=True
)

# lists containing mutations
medis = escaping_mutations.index.values
S309 = np.array(
    [escaping_mutations["Mutation"]["S309"]]
)  # ["S:I332","S:T333","S:N334","S:L335","S:C336","S:P337","S:G339","S:E340","S:V341","S:N343","S:A344","S:T345","S:R346","S:N354","S:K356","S:R357","S:I358","S:S359","S:N360","S:C361","S:N440","S:L441","S:R509"]
AZD1061 = np.array(
    [escaping_mutations["Mutation"]["AZD1061"]]
)  # ["S:T345","S:R346","S:N439","S:N440","S:L441","S:S443","S:K444","S:V445","S:G446","S:G447","S:Y449","S:N450","S:E484","S:F490","S:L492","S:Q493","S:S494","S:P499"]
AZD8895 = np.array(
    [escaping_mutations["Mutation"]["AZD8895"]]
)  # ["S:L455","S:F456","S:A475","S:G476","S:S477","S:T478","S:P479","S:E484","S:G485","S:F486","S:N487","S:C488","S:Y489","S:Q493"]

# create data structures for loop
df_sid = pd.DataFrame(
    columns=["Sample_id", "Antibody", "hig_imp_fac", "Escaping_mutations"]
)
df_g = pd.DataFrame(
    columns=["Sample_id", "Antibody", "hig_imp_fac", "Escaping_mutations"]
)
shortsigs = []

# creation of output file
for idx, (sid, group) in enumerate(df.groupby("Sample_id")):
    # check if there are any antibodies against which there are no mutations? <- Check with Folker
    group = group.reset_index(drop=True)
    for ix, sig in enumerate(group["Signature"]):
        short_sig = re.split("(\d+)", sig)[0] + re.split("(\d+)", sig)[1]
        shortsigs.append(short_sig)
    i_fac = factors[
        factors.index.str.contains("|".join(shortsigs), case=False, na=False)
    ]
    n_fac = i_fac.applymap(lambda x: (x - (-1000)) / (8.6 - (-1000)))
    ab = n_fac.sum().sort_values(ascending=False)
    for i, a in enumerate(ab.index):
        ls = sorted(
            list(zip(i_fac[ab.index[i]].index, i_fac[ab.index[i]])),
            key=lambda a: a[1],
            reverse=True,
        )
        hif = i_fac[a].min()
        df_g.at[i, "Escaping_mutations"] = ls
        df_g.at[i, "Antibody"] = a
        df_g.at[i, "hig_imp_fac"] = hif
        df_g.at[i, "Sample_id"] = sid
    df_sid = pd.concat([df_sid, df_g])

# for every sampleid check if present
# if present: fine
# else add with empty columns
df_sid = df_sid.reset_index(drop=True)
df_sid.to_json(snakemake.output[0], orient="records")
