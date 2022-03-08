from re import M

import pandas as pd
import pysam

AA_ALPHABET_TRANSLATION = {
    "Gly": "G",
    "Ala": "A",
    "Leu": "L",
    "Met": "M",
    "Phe": "F",
    "Trp": "W",
    "Lys": "K",
    "Gln": "Q",
    "Glu": "E",
    "Ser": "S",
    "Pro": "P",
    "Val": "V",
    "Ile": "I",
    "Cys": "C",
    "Tyr": "Y",
    "His": "H",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Thr": "T",
}


def df_from_vcf(vcf_file, variants=pd.DataFrame(columns=["Position", "Variant"])):

    with pysam.VariantFile(vcf_file, "r") as infile:
        for record in infile:
            for ann in record.info["ANN"]:
                ann = ann.split("|")
                hgvsp = ann[11]
                feature = ann[3]
                if hgvsp:
                    enssast_id, alteration = hgvsp.split(":", 1)
                    _prefix, alteration = alteration.split(".", 1)
                    for triplet, amino in AA_ALPHABET_TRANSLATION.items():
                        alteration = alteration.replace(triplet, amino)

                    hgvsp = f"{feature}:{alteration}"
                    # if alteration in snakemake.params.voc.get(feature, {}):
                    variants = variants.append(
                        {
                            "Position": str(record.pos),
                            "Variant": hgvsp,
                        },
                        ignore_index=True,
                    )
    # variants = variants.set_index("Variant")
    return variants


sanger_x_genome = pd.concat(
    [pd.read_csv(x) for x in snakemake.input.sanger_vs_ngs_genome]
)

print(sanger_x_genome)
min_pos = sanger_x_genome["Aln Start(t)"].min()
max_pos = sanger_x_genome["Aln End(t)"].max()
print(min_pos, max_pos)
sanger_x_genome.to_csv(snakemake.output.sanger_vs_genome)

sanger_variants = pd.DataFrame(columns=["Position", "Variant"])
for file in snakemake.input.sanger:
    sanger_variants = df_from_vcf(file, sanger_variants)

NGS_variants = df_from_vcf(snakemake.input.ngs_genome)

with open(snakemake.input.coverage, "r") as coverage_file:
    coverage = coverage_file.read().splitlines()

NGS_index_list = NGS_variants.set_index("Variant").index.tolist()
sanger_index_list = sanger_variants.set_index("Variant").index.tolist()

column = []
sanger_in_ngs = 0
for var in sanger_index_list:
    # print(var)
    # print(NGS_index_list)
    if var in NGS_index_list:
        sanger_in_ngs += 1
        column.append(True)
    else:
        column.append(False)

sanger_variants["in_ngs"] = column

with open(snakemake.output.sanger_vars, "w") as sanger_out:
    for index, var in sanger_variants.iterrows():
        print(index)
        print(coverage[int(var[0])])
        print(
            snakemake.wildcards.sample
            + ","
            + var[0]
            + ","
            + var[1]
            + ","
            + str(column[index])
            + ","
            + coverage[int(var[0])].split("\t")[-1],
            file=sanger_out,
        )

with open(snakemake.output.ngs_vars, "w") as ngs_out:
    for index, var in NGS_variants.iterrows():
        if int(var[0]) in range(int(min_pos), int(max_pos)):
            print(
                snakemake.wildcards.sample + "," + var[0] + "," + var[1], file=ngs_out
            )
print(NGS_variants)
column = []


with open(snakemake.output.sanger_ngs_diff, "w") as outfile, open(
    snakemake.output.sanger_ngs_diff_readable, "w"
) as outfile2:
    print(str(sanger_in_ngs) + "," + str(len(sanger_index_list)))
    print(str(sanger_in_ngs) + "," + str(len(sanger_index_list)), file=outfile)
    print(
        snakemake.wildcards.sample
        + ","
        + str(sanger_in_ngs)
        + "/"
        + str(len(sanger_index_list)),
        file=outfile2,
    )
