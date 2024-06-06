# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

# %%
sys.stderr = open(snakemake.log[0], "w")
import pysam

# see https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CINS = 1
BAM_CDEL = 2

N_WINDOW = 100


def is_uncertain_region(record, rpos, rend, ref_fasta):
    refseq = ref_fasta.fetch(record.reference_name, rpos - N_WINDOW, rend + N_WINDOW)
    return ("N" * 10) in refseq


def get_edit_dist(record, ref_fasta):
    edit_dist = 0
    qpos = 0
    rpos = record.reference_start
    for op, length in record.cigartuples:
        if op == BAM_CEQUAL:
            # no edit
            qpos += length
            rpos += length
        else:
            if op == BAM_CDIFF:
                refseq = ref_fasta.fetch(record.reference_name, rpos, rpos + length)
                if not is_uncertain_region(record, rpos, rpos + length, ref_fasta):
                    # Only consider edits where the reference has a true nucleotide
                    # because IUPAC codes lead to mason simulating an N in the reads.
                    edit_dist += sum(base in "ACGT" for base in refseq)
                qpos += length
                rpos += length

            elif op == BAM_CSOFT_CLIP:
                edit_dist += length
                qpos += length
            elif op == BAM_CHARD_CLIP:
                edit_dist += length
            elif op == BAM_CINS:
                refseq = ref_fasta.fetch(record.reference_name, rpos, rpos + length)
                if not is_uncertain_region(record, rpos, rpos + length, ref_fasta):
                    # only count edit distance if the region behind the insert does not
                    # contain N stretches in the reference. Rationale: those regions are apparently
                    # hard to assemble, and we cannot properly simulate reads for them, so
                    # we should not count them in a benchmark.
                    edit_dist += length
                qpos += length
            elif op == BAM_CDEL:
                refseq = ref_fasta.fetch(record.reference_name, rpos, rpos + length)
                if not is_uncertain_region(record, rpos, rpos + length, ref_fasta):
                    # only count edit distance if the region behind the insert does not
                    # contain N stretches in the reference. Rationale: those regions are apparently
                    # hard to assemble, and we cannot properly simulate reads for them, so
                    # we should not count them in a benchmark.
                    edit_dist += length
                rpos += length
            else:
                raise ValueError(f"Unsupported CIGAR operation: {op}")
    return edit_dist


with open(snakemake.output[0], "w") as out:
    print(
        "Accession",
        "Contigs",
        "Total contigs",
        "Contig length",
        "Reference length",
        "Contig frac",
        "Cum frac",
        "Relevant edit dist",
        "Cigar string",
        "Edit frac",
        sep="\t",
        file=out,
    )
    sum_of_edit_dist = 0
    largest_no_contig = 0
    small_larg_cover = 1.0
    small_larg_cover_name = "asd"
    smallest_cum_frac = 1.0

    # for bam_file in snakemake_input:
    for bam_file, ref_fasta in zip(snakemake.input.bams, snakemake.input.refs):
        current_contig = 1

        ref_fasta = pysam.FastaFile(ref_fasta)

        with pysam.AlignmentFile(bam_file, "rb") as samfile:
            total_contigs = samfile.count()
            accession = samfile.get_reference_name(0)

            largest_no_contig = (
                total_contigs
                if total_contigs > largest_no_contig
                else largest_no_contig
            )

        with pysam.AlignmentFile(bam_file, "rb") as samfile:
            ref_lengths = samfile.lengths[0]
            largest_coverage = 0
            cum_frac = 0
            for read in samfile.fetch():
                query_alignment_length = read.query_alignment_length
                frac = round(query_alignment_length / ref_lengths, 2)

                # Calculate edit distance from CIGAR string, because NM tag counts matching Ns as edits.
                edit = get_edit_dist(read, ref_fasta)

                cum_frac += frac
                cum_frac = round(cum_frac, 2)
                sum_of_edit_dist = sum_of_edit_dist + edit
                edit_frac = round(edit / query_alignment_length, 5)
                largest_coverage = frac if frac > largest_coverage else largest_coverage

                print(
                    accession,
                    current_contig,
                    total_contigs,
                    query_alignment_length,
                    ref_lengths,
                    frac,
                    cum_frac,
                    edit,
                    read.cigarstring,
                    edit_frac,
                    sep="\t",
                    file=out,
                )

                current_contig += 1
            smallest_cum_frac = (
                cum_frac if smallest_cum_frac > cum_frac else smallest_cum_frac
            )
            small_larg_cover = (
                largest_coverage
                if small_larg_cover > largest_coverage
                else small_larg_cover
            )
            small_larg_cover_name = (
                accession
                if small_larg_cover >= largest_coverage
                else small_larg_cover_name
            )

    print("Largest number of contigs")
    print(largest_no_contig)
    print("Smallest largest coverage in " + str(small_larg_cover_name))
    print(small_larg_cover)
    print("Sum of edit distances")
    print(sum_of_edit_dist)
    print("Smallest cum. fraction of contigs")
    print(smallest_cum_frac)
    print("Largest number of contigs", file=out)
    print(largest_no_contig, file=out)
    print("Smallest largest coverage in " + str(small_larg_cover_name), file=out)
    print(small_larg_cover, file=out)
    print("Sum of edit distances", file=out)
    print(sum_of_edit_dist, file=out)
    print("Smallest cum. fraction of contigs", file=out)
    print(smallest_cum_frac, file=out)
