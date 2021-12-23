import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from pysam import VariantFile

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


def phred_to_prob(phred):
    if phred is None:
        return 0
    return 10 ** (-phred / 10)


def get_AA_variant(record):
    for ann in record.info["ANN"]:
        ann = ann.split("|")
        hgvsp = ann[11]
        feature = ann[3]
        if hgvsp:
            _enssast_id, alteration = hgvsp.split(":", 1)
            _prefix, alteration = alteration.split(".", 1)
            for triplet, amino in AA_ALPHABET_TRANSLATION.items():
                alteration = alteration.replace(triplet, amino)
            hgvsp = f"{feature}:{alteration}"
    return hgvsp


def extract_data(variant_file: VariantFile):
    variants = []
    for record in variant_file.fetch():
        prob_not_present = phred_to_prob(record.info["PROB_ABSENT"][0]) + phred_to_prob(
            record.info["PROB_ARTIFACT"][0]
        )
        variant = get_AA_variant(record)
        if variant != "":
            variants.append(
                {
                    "chrom": record.chrom,
                    "pos": record.pos,
                    "variant": variant,
                    # "ANN" : record.info["ANN"],
                    "vaf": record.samples[0]["AF"][0],
                    "prob_not_present": prob_not_present,
                }
            )
    data = pd.DataFrame(variants)
    data["prob_present"] = 1 - data["prob_not_present"]
    data.drop(columns="prob_not_present", inplace=True)
    data.drop_duplicates(subset=["chrom", "pos", "variant"], inplace=True)
    return data


found_variants_list = []
for illumina_bcf, ont_bcf in zip(snakemake.input.illumina_bcf, snakemake.input.ont_bcf):
    with VariantFile(illumina_bcf, "rb") as illumina_file:
        illumina_variants = extract_data(illumina_file)

    with VariantFile(ont_bcf, "rb") as ont_file:
        ont_variants = extract_data(ont_file)

    found_variants = pd.merge(
        illumina_variants,
        ont_variants,
        how="outer",
        left_on=["chrom", "pos", "variant"],
        right_on=["chrom", "pos", "variant"],
        suffixes=("_illumina", "_ont"),
    )
    found_variants_list.append(found_variants)

found_variants = pd.concat(found_variants_list)

# found_variants.fillna({'vaf_illumina':0.0, 'vaf_ont':0.0, 'prob_present_illumina':0.0, 'prob_present_ont':0.0}, inplace=True)
found_variants["difference"] = found_variants["vaf_illumina"].astype(
    float
) - found_variants["vaf_ont"].astype(float)
found_variants.sort_values(by="difference", inplace=True, ascending=False)

found_variants.drop(columns=["difference"], inplace=True)
found_variants["test_case"] = snakemake.wildcards.test_case
found_variants.to_csv(snakemake.output[0], index=False, sep="\t")
