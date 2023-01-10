# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

# import required packages
import pathlib as path
import os
from xml.dom import INDEX_SIZE_ERR
import pandas as pd
import pysam
import re
import numpy as np
import requests
import glob
import sys

sys.stderr = open(snakemake.log[0], "w")

# transfrom phred score to probability
def phred_to_prob(phred):
    if phred is None:
        return pd.NA
    return 10 ** (-phred / 10)

# amino acid alphabet to translate 3 to 1 letter code
AA_ALPHABET_TRANSLATION = {
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

# create empty lists
data_dict = {
    "Signature": [],
    "Gene": [],
    "Position": [],
    "ReadDepth": [],
    "Frequency": [], 
    "Probability": [],
    "Sample_id": [],
    "Lineage":[]
}

# iterate over files in directory
for file in snakemake.input.bcfs:
    bcf_in = (pysam.VariantFile(file))
    lineages = pd.read_csv(snakemake.input.pan_calls, dtype={'Sample': str})
    idx = list(lineages['Sample']).index(str(filename))
    # pseudo or scaffold(preferred) from pangolin_calls_per_stage.csv if scaff empty take pseudo
    lineage = lineages['Scaffolded Seq'].values[idx] if lineages['Scaffolded Seq'].values[idx] != "Not called on" else lineages['Pseudo Seq'].values[idx]
    # fill lists with data obtained from file
    # iterate over records
    for record in bcf_in:
        pattern_count = str(record).count(':p.')
        # iterate over annotations in record
        # as there are multiple annotations per line, which cause problems, check for pattern (Aminoacid/protein annotation) first
        for ann in record.info["ANN"]:
            # split string seperated by pipe
            anno = ann.split("|")
            # get location
            loc = anno[3]
            # get feature_id
            fea_id = anno[6]
            # get alteration
            alt = anno[11]
            # check for occurence of pattern
            if pattern_count == 0 and alt != "":
                continue
            elif pattern_count != 0 and alt == "":
                continue
            elif pattern_count != 0 and alt != "":
                regex=re.compile(r':p.')
                alteration = regex.split(alt)[1]
                # replace %3D and triplet aa code with one letter code
                for triplet, amino in AA_ALPHABET_TRANSLATION.items():
                    alteration = alteration.replace("%3D", triplet[0]).replace(triplet, amino)
                alt = f"{loc}:{alteration}"
                pos = (record.pos)
                gene_name = loc
            elif pattern_count == 0 and alt == "":
                # if regular name is not available annotate using nucleotides
                ref=(record.ref)
                pos=(record.pos)
                # if clause to find mutation on ORFs
                if 266 <= pos <= 13468:
                    gene_name="ORF1a"
                elif 13468 <= pos <= 21555:
                    gene_name="ORF1b"
                elif 21563 <= pos <= 25384:
                    gene_name="S"
                elif 25393 <= pos <= 26220:
                    gene_name="ORF3a"
                elif 25765 <= pos <= 26220:
                    gene_name="ORF3b"
                elif 26245 <= pos <= 26472:
                    gene_name="E"
                elif 26523 <= pos <= 27191:
                    gene_name="M"
                elif 27202 <= pos <= 27387:
                    gene_name="ORF6"
                elif 27394 <= pos <= 27759:
                    gene_name="ORF7a"
                elif 27756 <= pos <= 27887:
                    gene_name="ORF7b"
                elif 27894 <= pos <= 28259:
                    gene_name="ORF8"
                elif 28284 <= pos <= 28577:
                    gene_name="ORF9a"
                elif 28734 <= pos <= 28955:
                    gene_name="ORF9b"
                elif 28274 <= pos <= 29533:
                    gene_name="N"
                elif 29558 <= pos <= 29674:
                    gene_name="ORF10"
                else:
                    gene_name="none"
                alts=(record.alts[0])
                alt=f"{ref}{pos}{alts}"
            else:
                raise Exception("Something unexpected occured in this record:", str(record), file)
            # add results in lists
            data_dict["Signature"].append(alt)
            # position
            data_dict["Position"].append(pos)
            # ORF
            data_dict["Gene"].append(gene_name)
            # read depth
            data_dict["ReadDepth"].append(record.samples[0]["DP"])
            # frequency
            data_dict["Frequency"].append(record.samples[0]["AF"][0])
            # probability
            data_dict["Probability"].append(1.0 - (phred_to_prob(record.info["PROB_ABSENT"][0]) + phred_to_prob(record.info["PROB_ARTIFACT"][0])))
            # filename
            data_dict["Sample_id"].append(filename)
            # lineage
            data_dict["Lineage"].append(lineage)    
    
# create dataframe
df = pd.DataFrame().from_dict(data_dict)

# filter duplicates 
df = df.drop_duplicates()

# export list of all mutations found
df.to_csv(snakemake.output, index=False, sep=",")
