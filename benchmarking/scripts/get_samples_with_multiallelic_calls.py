import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam

with open(snakemake.output[0], "w") as out_file:
    for vcf_path, meta in zip(snakemake.input.vcfs, snakemake.params.metadata):
        workflow, sample = meta.split(",")

        # get all refs and positions
        refs = []
        with pysam.VariantFile(vcf_path, "rb") as in_vcf:
            for record in in_vcf.fetch():
                refs.append(f"{record.ref},{record.pos}")

        # check if there are duplicated ref pos
        # these indicate multiallelic calls.
        if len(refs) != len(set(refs)):
            print(f"{sample}\t{workflow}\t{vcf_path}", file=out_file)
