import sys

sys.stderr = open(snakemake.log[0], "w")

import os
from datetime import date

min_length = int(snakemake.params.min_length)

length = 0
header2seq = {}
for file in snakemake.input.contigs:
    with open(file, "r") as infile:
        # TODO replace the following with a pysam one/two-liner
        for line in infile:
            line = line.rstrip()
            if line.startswith(">"):
                header2seq[">" + file.split("/")[-1].split(".")[0]] = ""
                length = 0
            else:
                header2seq[">" + file.split("/")[-1].split(".")[0]] += line
                length += len(line)

# TODO replace the following with either a pandas dataframe based code or with the
# Python csv module (delimiter=";")
outtab = open(snakemake.output.table, "a")
outfile = open(snakemake.output.fasta, "a")
if os.path.getsize(snakemake.output.table) == 0:
    outtab.write(
        "SENDING_LAB;DATE_DRAW;SEQ_TYPE;SEQ_REASON;SAMPLE_TYPE;PUBLICATION_STATUS;OWN_FASTA_ID\n"
    )
for key in header2seq:
    if len(header2seq[key]) > min_length:
        outtab.write("10259;;ILLUMINA;N;s001;N;%s\n" % (key[1:]))
        outfile.write(key + "\n" + header2seq[key] + "\n")
outfile.close()
outtab.close()
