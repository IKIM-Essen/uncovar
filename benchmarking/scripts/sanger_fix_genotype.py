import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam

with pysam.VariantFile(snakemake.input[0]) as in_vcf:
    print("--> input ", snakemake.input[0], file=sys.stderr)
    header = in_vcf.header

    with pysam.VariantFile(snakemake.output[0], "w", header=header) as out_vcf:
        # check each record for GT tag
        for record in in_vcf.fetch():
            for sample in record.samples.keys():
                print(f"Set GT tag of sample {sample} to 0/1", file=sys.stderr)
                record.samples[sample]["GT"] = (0, 1)

            out_vcf.write(record)
