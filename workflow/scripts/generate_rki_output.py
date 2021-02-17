import os
import sys
from datetime import date

min_length = int(snakemake.params.get("min_length"))

sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.log[0], "a")

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

outtab = open(snakemake.output.table, "a")
outfile = open(snakemake.output.fasta, "a")

if os.path.getsize(snakemake.output.table) == 0:
        outtab.write("IMS_ID;SENDING_LAB;DATE_DRAW;SEQ_TYPE;SEQ_REASON;SAMPLE_TYPE;OWN_FASTA_ID\n")
countlist = open("resources/countlist.txt", "r")
for line in countlist:
        counter = int(line.split(",")[1]) # get current start for numbering
countlist.close()
countlist = open("resources/countlist.txt", "a")
for key in header2seq:
        print(key)
        if len(header2seq[key]) > min_length:
                counter += 1
                outtab.write("IMS-10259-CVDP-%05d;IMS-10259-CVDP-%05d;;ILLUMINA;N;s001;%s\n" % (counter, counter, key[1:]))
                outfile.write(key + "\n" + header2seq[key] + "\n")
                countlist.write(key[1:] + "," + str(counter) + "\n")
                
outfile.close()
outtab.close()