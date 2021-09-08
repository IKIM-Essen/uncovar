sys.stderr = open(snakemake.log[0], "w")

import pysam


def extract_coverage_and_mask(
    bamfile_path: str,
    sequence_path: str,
    masked_sequence_path: str,
    coverage_path: str,
    min_coverage: int,
    min_allele: float,
    coverage_header: str = "#CHROM\tPOS\tCoverage\n",
) -> None:
    """Masks positions below a certain coverage with "N". Outputs the coverage per position in a separate file.

    Args:
        bamfile_path (str): Path to bamfile (.bam)
        sequence_path (str): Path to sequence (.fasta)
        sequence_path_index (str): Path to sequence index (.fasta.idx)
        masked_sequence_path (str): Path to write masked sequence to (.fasta)
        coverage_path (str): Path to write coverage per positions to (.txt)
        min_coverage (int): Minimal coverage at position to achive. Else apply masking with "N"
        min_allele (float): Minimal Allele frequency of base to achive. If not reach, masking base with "N".
        coverage_header (str, optional): Content of the header in the coverage file. Defaults to "#CHROM\tPOS\tCoverage\n".

    Raises:
        ValueError: if sequence contains more than one reference / contig.
    """
    # source: https://www.bioinformatics.org/sms/iupac.html
    IUPAC = {
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
        "B": ["C", "G", "T"],
        "D": ["A", "G", "T"],
        "H": ["A", "C", "T"],
        "V": ["A", "C", "G"],
        "N": ["A", "C", "T", "G"],
    }

    # sort iupac for later matching
    for key, values in IUPAC.items():
        IUPAC[key] = sorted(values)

    # context managers for bamfile reader, sequence reader and coverage writer
    with pysam.AlignmentFile(bamfile_path, "rb") as bamfile, open(
        sequence_path
    ) as sequence_handle, open(coverage_path, "w") as coverage:

        # get sequence(s) in fasta file
        # TODO replace FASTA parsing with pysam code
        sequence_dict = {}
        for line in sequence_handle:
            line = line.strip()
            if line.startswith(">"):
                key = line
                sequence_dict[key] = ""
            else:
                sequence_dict[key] += line

        if len(sequence_dict.keys()) > 1:
            raise ValueError("Sequence contains more than one contig.")

        # convert sequence string to list of characters so that we can change characters
        sequence = list(list(sequence_dict.values())[0])

        # write header of coverage file
        if len(coverage_header) > 0:
            coverage.write(coverage_header)

        # pileup reads per position
        for pileupcolumn in bamfile.pileup(ignore_overlaps=False):

            # write name, pos and coverage to coverage file
            coverage.write(
                "%s\t%s\t%s\n"
                % (
                    pileupcolumn.reference_name,
                    pileupcolumn.reference_pos,
                    pileupcolumn.nsegments,
                )
            )

            # check if there is enough coverage at poisition
            if pileupcolumn.nsegments < min_coverage:
                # log the masking
                print(
                    "Coverage of base %s at pos. %s = %s. Masking with N."
                    % (
                        sequence[pileupcolumn.reference_pos],
                        pileupcolumn.reference_pos,
                        pileupcolumn.nsegments,
                    ),
                    file=sys.stderr,
                )

                # mask the position
                sequence[pileupcolumn.reference_pos] = "N"

            # check Allele frequency of base
            else:
                # get all base in read for pileup base
                pileupread_base_count = {}
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:

                        read_base = pileupread.alignment.query_sequence[
                            pileupread.query_position
                        ]

                        if read_base not in pileupread_base_count.keys():
                            pileupread_base_count[read_base] = 1
                        else:
                            pileupread_base_count[read_base] += 1

                # calc Allele frequency
                try:
                    base_allel_freq = (
                        pileupread_base_count[sequence[pileupcolumn.reference_pos]]
                        / pileupcolumn.nsegments
                    )
                except:
                    base_allel_freq = 0

                # if Allele frequency is lower than treshold, then
                if base_allel_freq < min_allele:

                    bases = sorted(list(pileupread_base_count.keys()))

                    # mask base with corresonding IUPAC code
                    if len(bases) > 1:
                        for key, values in IUPAC.items():
                            if values == bases:
                                iupac_mask = key
                                break

                    # mask the 1 base with "N"
                    else:
                        iupac_mask = "N"

                    # log when masking occures
                    print(
                        "Coverage of base %s at pos. %s = %s with Allel frequency = %s. Bases in reads: %s. Masking with %s."
                        % (
                            sequence[pileupcolumn.reference_pos],
                            pileupcolumn.reference_pos,
                            pileupcolumn.nsegments,
                            base_allel_freq,
                            pileupread_base_count,
                            iupac_mask,
                        ),
                        file=sys.stderr,
                    )

                    sequence[pileupcolumn.reference_pos] = iupac_mask

    # join list of characters to sequence
    sequence = "".join(sequence)
    # TODO replace this mess with more clearer code
    header = list(sequence_dict.keys())[0].split(".")[0] + "\n"

    # write masked fasta file
    with open(masked_sequence_path, "w") as w:
        print(header, file=w)
        print(sequence, file=w)


if __name__ == "__main__":
    extract_coverage_and_mask(
        snakemake.input.bamfile,
        snakemake.input.sequence,
        snakemake.output.masked_sequence,
        snakemake.output.coverage,
        snakemake.params.min_coverage,
        snakemake.params.min_allele,
    )
