# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys

sys.stderr = open(snakemake.log[0], "w")
# sys.stdout = open(snakemake.log[0], "a")

import json

import pandas as pd
import pysam
import re

def iter_with_samples(inputfiles):
    return zip(snakemake.params.samples, inputfiles)

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

data = pd.DataFrame()

sample = snakemake.wildcards.sample
file = snakemake.input.bcf
variants_of_interest = {}
other_variants = {}

def insert_entry(variants, hgvsp, vaf, nucl_pos, pos, number):
    if variants.get(hgvsp):
        prev_vaf = variants.get(hgvsp)
        if prev_vaf[0] is None or prev_vaf[0] < vaf:
            # Only insert if there was no entry before or it had a smaller vaf.
            # Such duplicate calls can occur if there are multiple genomic variants
            # that lead to the same protein alteration.
            # We just report the protein alteration here, so what matters to us is the
            # variant call with the highest VAF.
            # TODO: in principle, the different alterations could even be complementary.
            # Hence, one could try to determine that and provide a joint vaf.
            variants[hgvsp] = [vaf, nucl_pos, pos, number]
    else:
        variants[hgvsp] = [vaf, nucl_pos, pos, number]

def fmt_variants(variants):
    return " ".join(sorted(f"{entry[1]} {hgvsp}:{entry[0]:.3f}" for hgvsp, entry in variants.items()))

with pysam.VariantFile(file, "rb") as infile:
    for record in infile:
        vaf = record.samples[0]["AF"][0]
        number = record.samples[0]["DP"][0]
        pos = record.pos
        for ann in record.info["ANN"]:
            ann = ann.split("|")
            hgvsp = ann[11]
            enssast_id = ann[6]
            feature = ann[3]
            nucl_pos = ann[10]
            bef_after = ann[16]
            if hgvsp:
                # TODO think about regex instead of splitting
                enssast_id, alteration = hgvsp.split(":", 1)
                _prefix, alteration = alteration.split(".", 1)
                for triplet, amino in AA_ALPHABET_TRANSLATION.items():
                    alteration = alteration.replace(triplet, amino)
                enssast_id, nuc_alteration = nucl_pos.split(":", 1)
                _prefix, nuc_alteration = nuc_alteration.split(".", 1)
                hgvsp = f"{feature}:{alteration}"
                print(nuc_alteration)
                
                if not "del" in nuc_alteration and not "ins" in nuc_alteration and not "inv" in nuc_alteration:
                    # reformat nucleotide position
                    m = re.match('([0-9]{1,})([ACTG])[>]([ACTG])', nuc_alteration)
                    nucl_pos = m.group(2) + m.group(1) + m.group(3)
                    insert_entry(variants_of_interest, hgvsp, vaf, nucl_pos, pos, number)
                elif "del" in nuc_alteration and not "ins" in nuc_alteration and not "inv" in nuc_alteration and "del" in hgvsp and not "ins" in hgvsp and not "inv" in hgvsp:
                    # reformat nucleotide position
                    coordinates = re.findall('([0-9]{1,})', nuc_alteration)
                    nucl_pos = "Deletion%s" % (str("-".join(coordinates))) 
                    # splitup multi AA deletions
                    feature, mutation = hgvsp.split(":", 1)
                    m = re.findall('([A-Z]{1}[0-9]{1,})', mutation)
                    for match in m:
                        hgvsp = feature + ":" + match + "-"
                        insert_entry(variants_of_interest, hgvsp, vaf, nucl_pos, pos, number)
                elif "delins" in nuc_alteration and not "inv" in nuc_alteration and not "inv" in hgvsp:
                    coordinates = re.findall('([0-9]{1,})', nuc_alteration)
                    insertion = re.search('[A-Z]+$', nuc_alteration)
                    nucl_pos = "DelIns%s%s" % (str("-".join(coordinates)), insertion.group(0))
                    # splitup multi AA deletions
                    feature, mutation = hgvsp.split(":", 1)
                    if "_" in mutation:
                        m = re.findall('([A-Z]{1}[0-9]{1,})', mutation)
                        insertion = re.search('[A-Z]+$', mutation)
                        for i in range(len(m)):
                            hgvsp = feature + ":" + m[i] + str(insertion)[i]
                            insert_entry(variants_of_interest, hgvsp, vaf, nucl_pos, pos, number)
                    else:
                        insert_entry(variants_of_interest, hgvsp, vaf, nucl_pos, pos, number)
                else:
                    nucl_pos = nuc_alteration
                    insert_entry(variants_of_interest, hgvsp, vaf, nucl_pos, pos, number)
                

                # insert_entry(variants_of_interest, hgvsp, vaf, nucl_pos)

for variant in variants_of_interest:
    data.loc[variant, "frequency"] = variants_of_interest[variant][0]
    data.loc[variant, "nucleotides"] = variants_of_interest[variant][1]
    data.loc[variant, "position"] = variants_of_interest[variant][2]
    data.loc[variant, "num_reads"] = variants_of_interest[variant][3]

with open("resources/voc-variants/Omicron", "r") as infile:
    for var in infile.read().splitlines():
        data.loc[var, "omicron"] = var

# data = data.sort_values(by=["position"])
data.sort_index(inplace=True)
data.to_csv(snakemake.output.var_table)
    