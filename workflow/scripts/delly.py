from snakemake.shell import shell
import pysam

exclude = (
    "-x {}".format(snakemake.input.exclude)
    if snakemake.input.get("exclude", "")
    else ""
)

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

with pysam.AlignmentFile(snakemake.input.sample) as bam:
    read_count = sum(1 for read in bam)
    if read_count < 10000:
        # Not enough reads to perform SV calling.
        # Output empty BCF.
        with pysam.VariantFile(snakemake.output[0], "wb", header=pysam.VariantHeader()) as bcf:
            pass

shell(
    "OMP_NUM_THREADS={snakemake.threads} delly call {extra} "
    "{exclude} -g {snakemake.input.ref} "
    "-o {snakemake.output[0]} {snakemake.input.sample} {log}"
)