from numpy import cov
import pandas as pd

for vars, covs in zip(snakemake.input.variants, snakemake.input.coverage):
    print(vars, covs)
    var_df = pd.read_csv(vars, header=None, index_col=1)
    cov_df = pd.read_csv(covs, sep="\t")

    print(var_df)
    print(cov_df)
        

with open(snakemake.output.table, "w") as outfile:
    print(" ", file=outfile)