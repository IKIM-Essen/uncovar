import pandas as pd
import sys
import pprint
import os
import pysam

sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.log[0], "a")
KRAKEN_FILTER_KRITERIA = "D"


final_df = pd.DataFrame()
for file in snakemake.input.contigs:
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
    final_df = final_df.append(
        {
            "contig length": str(length),
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
    krake_df_filtered["thereof SARS"] = krake_df[
        krake_df["name"] == "Severe acute respiratory syndrome-related coronavirus"
    ]["%"].values[0]
    krake_df_filtered["Unclassified"] = krake_df[
        krake_df["name"] == "unclassified"
    ]["%"].values[0]

    krake_df_filtered = krake_df_filtered.set_index("sample")
    total_kraken_df = total_kraken_df.append(krake_df_filtered)

# assembly_df = assembly_df.merge(
#     final_df, how="left", left_on="sample", right_on="sample"
# )

table = {}
for file in snakemake.input.pangolin:

    pang_call = open(file, "r")
    table[file.split("/")[-1].split(".")[0]] = [[],[],[]]
    for line in pang_call.read().splitlines():
        if not line.startswith("taxon"):
            table[file.split("/")[-1].split(".")[0]][0] = [line.split(",")[1] + " " + line.split(",")[-1].split(" ")[0]]

AS3to1 = {
    "Gly": "G",
    "Ala": "A",
    "Leu": "L",
    "Met": "M",
    "Phe": "F",
    "Trp": "W",
    "Lys": "K",
    "Gln": "Q",
    "Glu": "E",
    "Ser": "S",
    "Pro": "P",
    "Val": "V",
    "Ile": "I",
    "Cys": "C",
    "Tyr": "Y",
    "His": "H",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Thr": "T",
}

for file in snakemake.input.bcf:
    variants = pysam.VariantFile(file, "rb")
    for record in variants:
        vaf = record.samples[0]["AF"]
         #print(vaf)
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
                if  feature == "S" and\
                    ("501" in alt or \
                    ("484" in alt and not "6484" in alt ) or\
                    "417" in alt or\
                    "6970" in alt):
                    hgvsp = f"{feature}:{alt}"
                    table[file.split("/")[-1].split(".")[0]][1].append(hgvsp + " " + str(round(vaf[0], 3)) + " " + enssast_id)
                elif vaf[0] > 0.5:
                    table[file.split("/")[-1].split(".")[0]][2].append(hgvsp + " " + str(round(vaf[0], 3)))

var_df = pd.DataFrame()
for sample in table: 
    var_df = var_df.append(
        {
            "other variants": " ".join(table[sample][2]),
            "variants of interest": " ".join(table[sample][1]),
            "pangolin strain": " ".join(table[sample][0]),
            "sample": sample,
        },
        ignore_index=True,
        
    )
var_df = var_df.set_index("sample")
print(var_df)
var_df = var_df[["pangolin strain", "variants of interest", "other variants"]]


final_df = final_df.merge(
    total_kraken_df, how="left", left_on="sample", right_on=total_kraken_df.index
).rename(columns={"key_0": "sample"})

final_df = final_df.merge(
    var_df, how="right", left_on="sample", right_on=var_df.index
)

final_df = final_df.set_index("sample")

print(final_df)
final_df.to_csv(snakemake.output[0])