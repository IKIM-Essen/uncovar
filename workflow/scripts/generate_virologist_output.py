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
for file in snakemake.input.reads_unfiltered:
    
    sample = file.split("/")[-1].split(".")[0]
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
for file in snakemake.input.reads_filtered:
    infile = open(file, "r")
    sample = file.split("/")[-2]
    for line in infile:
        if "max length" in line:
            num_reads = line.split(", ")[-2].split(" ")[0]
            filtered_reads_df = filtered_reads_df.append(
                {
                    "# filtered reads": f"{int(num_reads):,}",
                    "sample": sample,
                },
                ignore_index=True,
            )
    infile.close()
filtered_reads_df = filtered_reads_df.set_index("sample")

initial_df = pd.DataFrame()
for file in snakemake.input.initial_contigs:
    contigs = {}
    sample = file.split("/")[-1].split(".")[0]
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
for file in snakemake.input.polished_contigs:
    contigs = {}
    sample = file.split("/")[-1].split(".")[0]
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

for file in snakemake.input.kraken:

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
    table[file.split("/")[-1].split(".")[0]] = [[],[],[],[],[],[],[],[],[],[],[],[]]
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
                if  feature == "S" and alt in snakemake.params.get("voc"):
                # if  feature == "S" and\
                    # ("501" in alt or \
                    # ("484" in alt and not "6484" in alt ) or\
                    # "417" in alt or\
                    # "6970" in alt):
                    table[file.split("/")[-1].split(".")[0]][1].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5 and feature == "S":
                    table[file.split("/")[-1].split(".")[0]][2].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5 and feature == "M":
                    table[file.split("/")[-1].split(".")[0]][3].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5 and feature == "E":
                    table[file.split("/")[-1].split(".")[0]][4].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5 and feature == "ORF1ab":
                    table[file.split("/")[-1].split(".")[0]][5].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5 and feature == "ORF3a":
                    table[file.split("/")[-1].split(".")[0]][6].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5 and feature == "ORF6":
                    table[file.split("/")[-1].split(".")[0]][7].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5 and feature == "ORF7a":
                    table[file.split("/")[-1].split(".")[0]][8].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5 and feature == "ORF8":
                    table[file.split("/")[-1].split(".")[0]][9].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5 and feature == "ORF10":
                    table[file.split("/")[-1].split(".")[0]][10].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5:
                    table[file.split("/")[-1].split(".")[0]][11].append(hgvsp + " " + str(round(vaf[0], 3)))
                
for sample in table:
    for i in range(1, len(table[sample])):
        hashing = {}
        for j in range(len(table[sample][i])):
            var = table[sample][i][j].split(" ")[0]
            freq = float(table[sample][i][j].split(" ")[1])
            print(var, freq)
            if var in hashing and freq > float(hashing[var].split(" ")[1]):
                hashing[var] = table[sample][i][j]
            elif var not in hashing:
                hashing[var] = table[sample][i][j]
        table[sample][i] = []
        for key in hashing:
            table[sample][i].append(hashing[key])

print(table)

var_df = pd.DataFrame()
for sample in table: 
    var_df = var_df.append(
        {
            "other variants": " ".join(table[sample][11]),
            "ORF10 variants": " ".join(table[sample][10]),
            "ORF8 variants": " ".join(table[sample][9]),
            "ORF7a variants": " ".join(table[sample][8]),
            "ORF6 variants": " ".join(table[sample][7]),
            "ORF3a variants": " ".join(table[sample][6]),
            "ORF1ab variants": " ".join(table[sample][5]),
            "E variants": " ".join(table[sample][4]),
            "M variants": " ".join(table[sample][3]),
            "other S variants": " ".join(table[sample][2]),
            "variants of interest": " ".join(table[sample][1]),
            "pangolin strain (#SNPs)": " ".join(table[sample][0]),
            "sample": sample,
        },
        ignore_index=True,
        
    )
var_df = var_df.set_index("sample")
var_df = var_df[["pangolin strain (#SNPs)", "variants of interest", "other S variants", "M variants", "E variants", "ORF1ab variants", "ORF3a variants", "ORF6 variants", "ORF7a variants", "ORF8 variants", "ORF10 variants", "other variants"]]

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
variant_df = output_df.copy()
qc_df.drop(columns=["other S variants", "M variants", "E variants", "ORF1ab variants", "ORF3a variants", "ORF6 variants", "ORF7a variants", "ORF8 variants", "ORF10 variants", "other variants"], inplace=True)

try:
    variant_df.drop(columns=["# raw reads", "# trimmed reads", "# filtered reads", "initial contig (bp)", "final contig (bp)", "Eukaryota (%)", "Bacteria (%)", "Viruses (%)", "thereof SARS (%)", "Unclassified (%)"], inplace=True)
except:
    variant_df.drop(columns=["# raw reads", "# trimmed reads", "# filtered reads", "initial contig (bp)", "final contig (bp)"], inplace=True)
output_df.to_csv(snakemake.output.all_data)
qc_df.to_csv(snakemake.output.qc_data)
variant_df.to_csv(snakemake.output.var_data)
