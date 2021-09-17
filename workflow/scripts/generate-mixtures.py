sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as out:
    print(*snakemake.params.mixtures, sep="\n", file=out)
