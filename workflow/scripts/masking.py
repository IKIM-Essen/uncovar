import sys

sys.stderr = open(snakemake.log[0], "w")

import pysam
from collections import Counter

# source: https://www.bioinformatics.org/sms/iupac.html
IUPAC = {
    frozenset("A"): "A",
    frozenset("C"): "C",
    frozenset("G"): "G",
    frozenset("T"): "T",
    frozenset("AG"): "R",
    frozenset("CT"): "Y",
    frozenset("GC"): "S",
    frozenset("AT"): "W",
    frozenset("GT"): "K",
    frozenset("AC"): "M",
    frozenset("CGT"): "B",
    frozenset("AGT"): "D",
    frozenset("ACT"): "H",
    frozenset("ACG"): "V",
    frozenset("ACTG"): "N",
}


def get_sequence():
    with pysam.FastxFile(snakemake.input.sequence) as fh:
        for entry in fh:
            return entry.sequence


def split(word):
    return [char for char in word]


def get_base_count(pileupcolumn):
    bases = []

    # pileupread: representation of a read aligned to a particular position in the reference sequence.
    for pileupread in pileupcolumn.pileups:
        # TODO Check pileupread for missing bases
        if not pileupread.is_del and not pileupread.is_refskip:

            read_base = pileupread.alignment.query_sequence[pileupread.query_position]

            bases.append(read_base)

    return Counter(bases)


def get_coverage(base_count):
    return sum(base_count.values())


def get_and_write_coverages_and_base_counts(
    coverage_header: str = "#CHROM\tPOS\tCoverage",
):
    with pysam.AlignmentFile(snakemake.input.bamfile, "rb") as bamfile, open(
        snakemake.output.coverage, "w"
    ) as coverage_manager:
        print(coverage_header, file=coverage_manager)
        coverages = {}
        base_counts = {}
        for base in bamfile.pileup():
            base_count = get_base_count(base)
            coverage = get_coverage(base_count)
            print(
                "%s\t%s\t%s"
                % (
                    base.reference_name,
                    base.reference_pos,
                    coverage,
                ),
                file=coverage_manager,
            )

            coverages[base.reference_pos] = coverage
            base_counts[base.reference_pos] = base_count

        return coverages, base_counts


def get_allel_freq(correct_base, base_counts):
    return base_counts[correct_base] / sum(base_counts.values())


def get_UPAC_mask(base_counts):
    return IUPAC[frozenset("".join(base_counts.keys()))]


def mask_sequence(sequence, coverages, base_counts):
    sequence = split(sequence)

    covered_postions = coverages.keys()

    for position, base in enumerate(sequence):

        if position not in covered_postions:
            # TODO Check why there are postions that are not covered by any reads and are not Ns
            # sequence[position] = "N"

            print(
                "Base %s at pos. %s not covered by any read."
                % (
                    base,
                    position,
                ),
                file=sys.stderr,
            )
            continue

        if coverages[position] < snakemake.params.min_coverage:
            sequence[position] = "N"

            print(
                "Coverage of base %s at pos. %s = %s. Masking with N."
                % (
                    base,
                    position,
                    coverages[position],
                ),
                file=sys.stderr,
            )
            continue

        if (
            not snakemake.params.is_ont
            and get_allel_freq(base, base_counts[position])
            < snakemake.params.min_allele
        ):
            if "N" in base_counts[position].keys():
                mask = "N"
            else:
                mask = get_UPAC_mask(base_counts[position])

            print(
                "Coverage of base %s at pos. %s = %s with Allel frequency = %s. Bases in reads: %s. Masking with %s."
                % (
                    base,
                    position,
                    coverages[position],
                    get_allel_freq(base, base_counts[position]),
                    base_counts[position],
                    mask,
                ),
                file=sys.stderr,
            )

            sequence[position] = mask

    return "".join(sequence)


def write_sequence(sequence):
    with pysam.FastxFile(snakemake.input.sequence) as infile, open(
        snakemake.output.masked_sequence, mode="w"
    ) as outfile:
        print(">%s" % next(infile).name.split(".")[0], file=outfile)
        print(sequence, file=outfile)


sequence = get_sequence()
assert isinstance(sequence, str), "More than one sequence in .fasta file."
coverages, base_counts = get_and_write_coverages_and_base_counts()
masked_sequence = mask_sequence(sequence, coverages, base_counts)
write_sequence(masked_sequence)
