import sys

sys.stderr = open(snakemake.log[0], "w")
#parameter = snakemake.params.get("parameter", "")

import random

def select_random_lineages(sm_input, sm_output, number_of_samples):

    with open (sm_input) as f:
        lines = f.read().splitlines()

    lineages = []

    for i in range(number_of_samples):
        rnd_strain_path = random.choice(lines)
        strain = rnd_strain_path.replace(".fasta", "").split("/")[-1]
        strain = strain.replace(".", "-")
        lineages.append(strain)

    with open(sm_output, "w") as handler:
        for element in lineages:
            handler.write(element + "\n")



if __name__ == "__main__":
    select_random_lineages(snakemake.input[0], snakemake.output[0], snakemake.params.number_of_samples)