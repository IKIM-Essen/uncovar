import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam

with pysam.VariantFile(snakemake.input[0]) as in_vcf, open(
    snakemake.output[0], "w"
) as out_file:
    header = in_vcf.header
    print(len([x for x in in_vcf.fetch()]), file=out_file)
