import re
import pysam
import numpy as np
import sys

sys.stderr = open(snakemake.log[0], "w")

IUPAC = {
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


def phred_to_prob(phred):
    if phred is None:
        return 0
    return 10 ** (-phred / 10)


with pysam.FastaFile(snakemake.input.fasta) as infasta, pysam.VariantFile(
    snakemake.input.bcf, "rb"
) as invcf, pysam.AlignmentFile(snakemake.input.bam, "rb") as inbam:

    assert len(infasta.references) == 1, "expected reference with single contig"
    contig = infasta.references[0]
    ref_seq = infasta.fetch(contig)
    cov_a, cov_c, cov_g, cov_t = inbam.count_coverage(contig)
    coverage = np.add(np.add(np.add(cov_a, cov_c), cov_g), cov_t)

    seq = ""
    last_pos = -1  # last considered reference position
    for record in invcf:
        rec_pos = record.pos - 1  # convert to zero based
        if rec_pos > last_pos + 2:
            chunk_seq = np.array(list(ref_seq[last_pos + 1 : rec_pos]))

            # check for low coverage regions
            chunk_low_cov = (
                coverage[last_pos + 1 : rec_pos] < snakemake.params.min_coverage
            )

            if len(chunk_seq[chunk_low_cov]) > 0:
                chunk_seq[chunk_low_cov] = "N"

            seq += "".join(chunk_seq)
        elif rec_pos < last_pos:
            # This must be an alternative allele to the last considered record.
            # But the last considered record had at least VAF>=0.5.
            # Hence, this must be a minor allele, and can therefore be ignored.
            continue

        last_pos = rec_pos - 1

        dp_sample = record.samples[0]["DP"][0]
        if dp_sample is None:
            dp_sample = 0

        # ignore low coverage records (subsequent iteration will add an N for that locus then)
        is_low_coverage = dp_sample < snakemake.params.min_coverage
        if is_low_coverage:
            continue

        def get_prob(event):
            return phred_to_prob(record.info[f"PROB_{event.upper()}"][0])

        prob_high = get_prob("clonal") + get_prob("subclonal_high")
        prob_major = get_prob("subclonal_major")

        apply = prob_high >= snakemake.params.min_prob_apply
        uncertain = (
            prob_major >= snakemake.params.min_prob_apply
            or (prob_high + prob_major) >= 0.5
        )

        if not (apply or uncertain):
            # we simply ignore this record
            continue

        assert len(record.alleles) == 2
        ref_allele, alt_allele = record.alleles

        # REF: A, ALT: <DEL>

        def handle_deletion(del_len):
            global last_pos
            global seq
            if not apply:
                seq += "N" * del_len
            last_pos += del_len + 1

        if alt_allele == "<DEL>":
            seq += ref_allele
            del_len = record.info["SVLEN"][0]
            handle_deletion(del_len)
        elif alt_allele == "<DUP>":
            dup_seq = ref_seq[rec_pos : record.stop]
            seq += dup_seq * 2
            last_pos += len(dup_seq)
        elif re.match("[A-Z]+$", alt_allele) is None:
            # TODO cover more variant types before publication
            raise ValueError(f"Unexpected alt allele: {alt_allele} not yet supported")
        elif len(ref_allele) == len(alt_allele):
            # SNV or MNV
            if apply:
                seq += alt_allele
            else:
                # store IUPAC codes
                for a, b in zip(*record.alleles):
                    bases = frozenset((a.upper(), b.upper()))
                    if len(bases) > 1:
                        # get IUPAC representation of bases
                        seq += IUPAC[bases]
                    else:
                        # add single base
                        (base,) = bases
                        seq += base
            last_pos += len(alt_allele)
        elif len(ref_allele) > 1 and len(alt_allele) == 1:
            # deletion
            del_len = len(ref_allele) - 1
            seq += alt_allele
            handle_deletion(del_len)
        elif len(ref_allele) == 1 and len(alt_allele) > 1:
            # insertion
            ins_seq = alt_allele[1:]
            seq += ref_allele
            if apply:
                seq += ins_seq
            else:
                seq += "N" * len(ins_seq)
            last_pos += 1
        elif len(ref_allele) > 1 and len(alt_allele) > 1:
            # replacement
            last_pos += len(ref_allele)
            if apply:
                seq += alt_allele
            else:
                seq += "N" * len(alt_allele)
        else:
            raise ValueError(f"Unexpected alleles: {ref_allele}, {alt_allele}")

    # add sequence until end
    seq += ref_seq[last_pos:]


with open(snakemake.output[0], "w") as outfasta:
    print(f">{snakemake.params.sample}", file=outfasta)
    print(seq, file=outfasta)
