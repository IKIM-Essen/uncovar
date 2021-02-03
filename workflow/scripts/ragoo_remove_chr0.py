from Bio import SeqIO
import sys
sys.stderr = open(snakemake.log[0], "w")

def remove_chr0(data_path, out_path):
    valid_records = []
    with open(data_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if 'Chr0' not in record.name:
                valid_records.append(record)
    SeqIO.write(valid_records, out_path, "fasta")

remove_chr0(snakemake.input[0], snakemake.output[0])