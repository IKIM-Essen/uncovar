import pysam
import os
import re

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

# variants = pysam.VariantFile(snakemake.input[0])
in_path = "/home/alex/repo/snakemake-workflow-sars-cov2/results/filtered-calls/ref~main/"
var_list = [f for f in os.listdir(in_path) if f.endswith('-impact.bcf')]
for file in var_list:
    # os.system("bcftools index " + in_path + file)
    variants = pysam.VariantFile(in_path + file, "rb")

    for record in variants:
        vaf = record.samples[0]["AF"]
         #print(vaf)
        for ann in record.info["ANN"]:
            ann = ann.split("|")
            hgvsp = ann[11]
            
            feature = ann[3]
            #print(feature)
            if hgvsp:
                alt = hgvsp.split(":", 1)[1].split(".")[1]
                match = re.match(r"([A-z]{3})([A-z0-9]+)([A-z]{3})", alt)
                #print(match)
                #alt = str(AS3to1[match[0]]) + str(match[1]) + str(AS3to1[match[2]])
                 # generate a set over all ann's, or, only take the one with the lexicographically lowest ENSSAST-ID (column 6) and put in table
                hgvsp = f"{feature}:{alt}"
                if  feature == "S" and\
                    ("501" in alt or \
                    ("484" in alt and not "6484" in alt ) or\
                    "417" in alt or\
                    "6970" in alt):
                    alt = str(AS3to1[match.group(1)] + str(match.group(2)) + str(AS3to1[match.group(3)]))
                    hgvsp = f"{feature}:{alt}"
                    print(file.split(".")[0],hgvsp, vaf)