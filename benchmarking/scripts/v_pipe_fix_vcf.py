import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam


def write_record(outvcf, new_record):
    record = outvcf.new_record()
    record.alleles = new_record.alleles
    record.chrom = new_record.chrom
    record.contig = new_record.contig
    record.pos = new_record.pos
    record.qual = new_record.qual
    record.id = new_record.id
    record.filter.add(new_record.filter.keys()[0])
    record.info.update(new_record.info)

    outvcf.write(record)


with pysam.VariantFile(snakemake.input[0]) as in_vcf:
    header = in_vcf.header
    header.info.add(
        id="SVLEN", number="1", type="Integer", description="Length of deletion"
    )
    with pysam.VariantFile(snakemake.output[0], "w", header=header) as out_vcf:
        for record in in_vcf.fetch():
            if record.alts == ("-",):
                record.alts = ("<DEL>",)
                record.info.update({"SVLEN": -len(record.alts)})

            write_record(out_vcf, record)
