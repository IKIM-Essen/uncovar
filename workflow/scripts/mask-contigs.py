sys.stderr = open(snakemake.log[0], "w")

import pysam


def extract_coverage_and_mask(
    bamfile_path: str,
    sequence_path: str,
    masked_sequence_path: str,
    coverage_path: str,
    min_coverage: int,
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
        coverage_header (str, optional): Content of the header in the coverage file. Defaults to "#CHROM\tPOS\tCoverage\n".

    Raises:
        ValueError: if sequence contains more than one reference / contig.
    """

    # context managers for bamfile reader, sequence reader and coverage writer
    with pysam.AlignmentFile(bamfile_path, "rb") as bamfile, open(
        sequence_path
    ) as sequence_handle, open(coverage_path, "w") as coverage:

        # get sequence(s)
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

        # convert sequence string to list of characters
        sequence = list(list(sequence_dict.values())[0])

        if len(coverage_header) > 0:
            coverage.write(coverage_header)

        # pileup reade per position
        for pileupcolumn in bamfile.pileup():

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
                # log when masking occures
                print(
                    "Masked at pos",
                    pileupcolumn.reference_pos,
                    "in",
                    pileupcolumn.reference_name,
                    file=sys.stderr,
                )

                # mask the position
                sequence[pileupcolumn.reference_pos] = "N"

    # join list of characters to sequence
    sequence = "".join(sequence)
    header = list(sequence_dict.keys())[0].split(".")[0] + "\n"

    # write masked fasta file
    with open(masked_sequence_path, "w") as w:
        w.write(header), w.write(sequence)


if __name__ == "__main__":
    extract_coverage_and_mask(
        snakemake.input.bamfile,
        snakemake.input.sequence,
        snakemake.output.masked_sequence,
        snakemake.output.coverage,
        snakemake.params.get("min_coverage", ""),
    )
