import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam

FILLER_SAMPLE_NAME = "SAMPLE"


def copy_record(old_record, new_record):
    new_record.alleles = old_record.alleles
    new_record.chrom = old_record.chrom
    new_record.contig = old_record.contig
    new_record.pos = old_record.pos
    new_record.qual = old_record.qual
    new_record.id = old_record.id
    new_record.info.update(old_record.info)

    return new_record


with pysam.VariantFile(snakemake.input[0]) as in_vcf:
    print("--> input ", snakemake.input[0], file=sys.stderr)
    header = in_vcf.header

    # add GT if not in header
    if "GT" not in header.formats.keys():
        print("Added GT to header...", file=sys.stderr)
        header.formats.add(id="GT", number="1", type="String", description="Genotype")

    # add sample if vcf does not contain one
    if header.samples.__len__() == 0:
        header.add_sample(FILLER_SAMPLE_NAME)

    with pysam.VariantFile(snakemake.output[0], "w", header=header) as out_vcf:
        # check each record for GT tag
        for record in in_vcf.fetch():
            if "GT" not in record.format.keys():

                # if the vcf does not contain a sample, the record also doesn't have on
                # thus we need to create a new one, which contains a sample
                if len(record.samples.items()) == 0:
                    new_record = out_vcf.new_record()
                    record = copy_record(record, new_record)

                for sample in record.samples.keys():
                    print(f"Add GT tag to sample {sample}.", file=sys.stderr)
                    record.samples[sample]["GT"] = (0, 1)

            out_vcf.write(record)
