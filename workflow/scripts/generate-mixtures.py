sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as out:
    for mix in snakemake.params.mixtures:
        out.write(f"{mix}\n")
