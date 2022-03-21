import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam

with open(snakemake.output[0], "w") as outfile:
    for vcf_path, sample in zip(snakemake.input, snakemake.params.samples):
        with pysam.VariantFile(vcf_path, "rb") as in_vcf:
            count_of_calls = 0
            for record in in_vcf.fetch():
                count_of_calls += 1

            if count_of_calls == 0:
                print(sample, file=outfile)
