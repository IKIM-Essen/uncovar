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

def df_from_vcf(vcf_file, variants = pd.DataFrame(columns=["Position", "Variant"])):
    
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

sanger_variants = pd.DataFrame(columns=["Position", "Variant"])
for file in snakemake.input.sanger:
    sanger_variants = df_from_vcf(file, sanger_variants)
print(sanger_variants)
sanger_index_list = sanger_variants.index.tolist()
NGS_variants = df_from_vcf(snakemake.input.genome)
NGS_index_list = NGS_variants.index.tolist()
print(NGS_variants)
column = []
for var in sanger_index_list:
    if var in NGS_index_list:
        column.append("both")
    else:
        column.append("sanger_only")
        print("SANGER_ONLY")

sanger_variants["compare"] = column
sanger_variants = sanger_variants.set_index("Position")
sanger_variants.to_csv(snakemake.output[0], mode="a")
