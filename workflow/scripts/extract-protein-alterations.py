import pysam
import os
import pprint

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
table = {}
# variants = pysam.VariantFile(snakemake.input[0])
in_path_vcf = "/home/alex/repo/snakemake-workflow-sars-cov2/results/2021-03-01/filtered-calls/ref~main/"
var_list = [f for f in os.listdir(in_path_vcf) if f.endswith('high+moderate-impact.bcf')]
in_path_pangolin = "/home/alex/repo/snakemake-workflow-sars-cov2/results/2021-03-01/tables/strain-calls/"
for file in var_list:

    variants = pysam.VariantFile(in_path_vcf + file, "rb")
    pang_call = open(in_path_pangolin + file.split(".")[0] + ".strains.pangolin.csv", "r")
    table[file.split(".")[0]] = [[],[],[]]
    for line in pang_call.read().splitlines():
        if not line.startswith("taxon"):
            #print(file.split(".")[0], line.split(",")[1], line.split(",")[-1].split(" ")[0], end="\t")
            table[file.split(".")[0]][0] = [line.split(",")[1] + " " + line.split(",")[-1].split(" ")[0]]
            #table[file.split(".")[0]][0] = " ".join([line.split(",")[1], line.split(",")[-1].split(" ")[0]])

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
                for triplet, amino in AS3to1.items():
                    alt = alt.replace(triplet, amino)
                hgvsp = f"{feature}:{alt}"
                if  feature == "S" and\
                    ("501" in alt or \
                    ("484" in alt and not "6484" in alt ) or\
                    "417" in alt or\
                    "6970" in alt):
                    hgvsp = f"{feature}:{alt}"
                    table[file.split(".")[0]][1].append(hgvsp + " " + str(round(vaf[0], 3)))
                elif vaf[0] > 0.5:
                    table[file.split(".")[0]][2].append(hgvsp + " " + str(round(vaf[0], 3))) 

pprint.pprint(table)