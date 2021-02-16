import os
import sys
from datetime import date

sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.log[0], "a")
# today = date.today()


# current_date = today.strftime("%Y-%m-%d")
# input = snakemake.input[0]

# all_final_contigs = os.listdir(in_path)
length = 0

header2seq = {}
for file in snakemake.input:
        print(file)
        infile = open(file, "r")
        for line in infile.read().splitlines():
                if line.startswith(">"):
                        header2seq[">" + file.split("/")[-1].split(".")[0]] = ""
                        length = 0
                else:
                        header2seq[">" + file.split("/")[-1].split(".")[0]] += line
                        length += len(line)
        infile.close()

#outfile = open(out_path + current_date + "_uk-essen_rki.fasta", "a")
# outtab_path = out_path + current_date + "_uk-essen_rki.csv"
outtab = open(snakemake.output.table[0], "a")
outfile = open(snakemake.output.fasta[0], "a")
if os.path.getsize(snakemake.output.table[0]) == 0:
        outtab.write("IMS_ID;SENDING_LAB;DATE_DRAW;SEQ_TYPE;SEQ_REASON;SAMPLE_TYPE;OWN_FASTA_ID\n")
counter = 21  # current start for numbering
for key in header2seq:
        print(key)
        if len(header2seq[key]) > 15:
                outtab.write("IMS-10259-CVDP-%05d;IMS-10259-CVDP-%05d;;ILLUMINA;N;s001;%s\n" % (counter, counter, key[1:]))
                outfile.write(key + "\n" + header2seq[key] + "\n")
                counter += 1
outfile.close()
outtab.close()