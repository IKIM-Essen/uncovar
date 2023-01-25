import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
from pysam import VariantFile


def read_vcf_file(path):
    if path.endswith(".bcf"):
        mode = "rb"
    elif path.endswith(".vcf"):
        mode = "r"
    else:
        raise TypeError(f"File extension of path {path} not recognized.")

    with VariantFile(path, mode) as variant_file:
        variants = []
        for record in variant_file.fetch():
            print(record.alleles)
            REF, ALT = record.alleles
            variants.append({"REF": REF, "ALT": ALT, "POS": record.pos})
    return pd.DataFrame(variants)


first_pipeline = read_vcf_file(snakemake.input.fist_pipeline)
second_pipeline = read_vcf_file(snakemake.input.second_pipeline)

print(first_pipeline)
print(second_pipeline)
