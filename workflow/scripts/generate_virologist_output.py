import pandas as pd
import sys
import os
import pysam
import json

sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.log[0], "a")
KRAKEN_FILTER_KRITERIA = "D"

print(snakemake.params.get("voc"))

initial_reads_df = pd.DataFrame()
for sample, file in zip(snakemake.params.samples, snakemake.input.reads_unfiltered):
    with open(file) as read_json:
        number_reads = json.load(read_json)
    raw_reads = int(number_reads["summary"]["before_filtering"]["total_reads"])
    trimmed_reads = int(number_reads["summary"]["after_filtering"]["total_reads"])
    initial_reads_df = initial_reads_df.append(
        {
            "# raw reads": f"{raw_reads:,}",
            "# trimmed reads": f"{trimmed_reads:,}",
            "sample": sample,
        },
        ignore_index=True,
    )
initial_reads_df = initial_reads_df.set_index("sample")
initial_reads_df = initial_reads_df[["# raw reads", "# trimmed reads"]]

filtered_reads_df = pd.DataFrame()
for sample, file in zip(snakemake.params.samples, snakemake.input.reads_used_for_assembly):
    infile = open(file, "r")
    for line in infile.read().splitlines():
        try:
            num_reads = int(line)/4*2
        except:
            num_reads = 0
        filtered_reads_df = filtered_reads_df.append(
            {
                "# used reads": f"{int(num_reads):,}",
                "sample": sample,
            },
            ignore_index=True,
        )
    infile.close()
filtered_reads_df = filtered_reads_df.set_index("sample")

initial_df = pd.DataFrame()
for sample, file in zip(snakemake.params.samples, snakemake.input.initial_contigs):
    contigs = {}
    if os.stat(file).st_size == 0:
        contigs[sample] = ""
    else:
        with open(file, "r") as fasta_unordered:
             
            for line in fasta_unordered.read().splitlines():
                if line.startswith(">"):
                    key = line
                    contigs[key] = ""
                else:  
                    contigs[key] += line

    length_initial = 0
    for key in contigs:
        if len(contigs[key]) > length_initial:
            length_initial = len(contigs[key])

    initial_df = initial_df.append(
        {
            "initial contig (bp)": f"{int(length_initial):,}",
            "sample": sample,
        },
        ignore_index=True,
    )

final_df = pd.DataFrame()
for sample, file in zip(snakemake.params.samples, snakemake.input.polished_contigs):
    contigs = {}
    if os.stat(file).st_size == 0:
        contigs[sample] = ""
    else:
        with open(file, "r") as fasta_ordered:
            
            for line in fasta_ordered.read().splitlines():
                if line.startswith(">"):
                    contigs[sample] = ""
                else:  
                    contigs[sample] += line

    for key in contigs:
        length = len(contigs[key])
    length = f'{length:,}'
    final_df = final_df.append(
        {
            "final contig (bp)": length,
            "sample": sample,
        },
        ignore_index=True,
    )

total_kraken_df = pd.DataFrame()

for sample, file in zip(snakemake.params.samples, snakemake.input.kraken):

    sample = file.split("/")[-1].split(".")[0]
    krake_df = pd.read_csv(
        file,
        delimiter="\t",
        header=None,
    )
    krake_df.columns = ["%", "covered", "assigned", "code", "taxonomic ID", "name"]
    krake_df.name = krake_df.name.str.strip()
    krake_df_filtered = (
        krake_df[(krake_df["code"] == KRAKEN_FILTER_KRITERIA)][["%", "name"]]
        .set_index("name")
        .T
    )
    krake_df_filtered["sample"] = sample
    try:
        krake_df_filtered["thereof SARS"] = krake_df[
            krake_df["name"] == "Severe acute respiratory syndrome-related coronavirus"
        ]["%"].values[0]
    except:
        krake_df_filtered["thereof SARS"] = 0
    krake_df_filtered["Unclassified"] = krake_df[
        krake_df["name"] == "unclassified"
    ]["%"].values[0]

    krake_df_filtered = krake_df_filtered.rename(columns={"Eukaryota": "Eukaryota (%)", "Bacteria": "Bacteria (%)", "Viruses": "Viruses (%)", "thereof SARS": "thereof SARS (%)", "Unclassified": "Unclassified (%)"})
    krake_df_filtered = krake_df_filtered.set_index("sample")
    total_kraken_df = total_kraken_df.append(krake_df_filtered)

try:
    total_kraken_df.drop(columns=["Archaea"], inplace=True)
except:
    pass

table = {}
for file in snakemake.input.pangolin:

    pang_call = open(file, "r")
    table[file.split("/")[-1].split(".")[0]] = [[] for _ in range(12)]
    for line in pang_call.read().splitlines():
        if not line.startswith("taxon"):
            if line.split(",")[1].startswith("None"):
                table[file.split("/")[-1].split(".")[0]][0] = ["no strain called"]
            else:
                table[file.split("/")[-1].split(".")[0]][0] = [line.split(",")[1] + " " + "(" + line.split(",")[-1].split(" ")[0]+ ")"]

AS3to1 = {
    "Gly": "G", "Ala": "A", "Leu": "L", "Met": "M",
    "Phe": "F", "Trp": "W", "Lys": "K", "Gln": "Q",
    "Glu": "E", "Ser": "S", "Pro": "P", "Val": "V",
    "Ile": "I", "Cys": "C", "Tyr": "Y", "His": "H",
    "Arg": "R", "Asn": "N", "Asp": "D", "Thr": "T",
}

for file in snakemake.input.bcf:
    variants = pysam.VariantFile(file, "rb")
    for record in variants:
        vaf = record.samples[0]["AF"]
        for ann in record.info["ANN"]:
            ann = ann.split("|")
            hgvsp = ann[11]
            enssast_id = ann[6]
            feature = ann[3]
            if hgvsp:
                alt = hgvsp.split(":", 1)[1].split(".")[1]
                enssast_id = hgvsp.split(":", 1)[0]
                for triplet, amino in AS3to1.items():
                    alt = alt.replace(triplet, amino)
                hgvsp = f"{feature}:{alt}"
                entry = f"{hgvsp}:{vaf[0]:.3f}"
                sample = file.split("/")[-1].split(".")[0]
                if feature == "S" and alt in snakemake.params.get("voc"):
                    table[sample][1].append(entry)
                else:
                    table[sample][2].append(entry)

for sample in table:
    for i in range(1, len(table[sample])):
        hashing = {}
        for j in range(len(table[sample][i])):
            print(table[sample][i][j])
            var = ':'.join(table[sample][i][j].split(":")[:2])
            freq = float(table[sample][i][j].split(":")[2])
            print(var, freq)
            if var not in hashing or freq > float(hashing[var].split(":")[2]):
                hashing[var] = table[sample][i][j]
        table[sample][i] = []
        for element in hashing.values():
            table[sample][i].append(element)

var_df = pd.DataFrame()
for sample in table: 
    var_df = var_df.append(
        {
            "other variants": " ".join(table[sample][2]),
            "variants of interest": " ".join(table[sample][1]),
            "pangolin strain (#SNPs)": " ".join(table[sample][0]),
            "sample": sample,
        },
        ignore_index=True,
        
    )
var_df = var_df.set_index("sample")
var_df = var_df[["pangolin strain (#SNPs)", "variants of interest", "other variants"]]

initial_df = initial_df.set_index("sample")
final_df = final_df.set_index("sample")

output_df = pd.merge(
    initial_reads_df, filtered_reads_df, how="left", left_on="sample", right_on=initial_reads_df.index
)
output_df = output_df.merge(
    initial_df, how="left", left_on="sample", right_on=initial_df.index
)

output_df = output_df.merge(
    final_df, how="left", left_on="sample", right_on=initial_df.index
)
try:
    total_kraken_df = total_kraken_df[["Eukaryota (%)", "Bacteria (%)", "Viruses (%)", "thereof SARS (%)", "Unclassified (%)"]]
except:
    pass
output_df = output_df.merge(
    total_kraken_df, how="left", left_on="sample", right_on=total_kraken_df.index
).rename(columns={"key_0": "sample"})

output_df = output_df.merge(
    var_df, how="right", left_on="sample", right_on=var_df.index
)

output_df = output_df.set_index("sample")

qc_df = output_df.copy()

output_df.to_csv(snakemake.output.all_data)
qc_df.to_csv(snakemake.output.qc_data)
