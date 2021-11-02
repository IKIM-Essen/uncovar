present_vars = 0
overall_vars = 0

for file in snakemake.input:
    with open(file) as infile:
        vars_in_both, of_vars_sang = infile.read().splitlines()[0].split(",")
        present_vars += int(vars_in_both)
        overall_vars += int(of_vars_sang)

with open(snakemake.output[0], "w") as outfile:
    print("Of", overall_vars, "variants of concern mutations,", present_vars, "are present in both =", present_vars/overall_vars*100, "%", file=outfile)
