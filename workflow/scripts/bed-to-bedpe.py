# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import pandas as pd
import urllib.request

# Function to check bed for url
def check_bed_for_URL(bed_file):
    if "https" in bed_file:
        filename = bed_file.split("/")[-1]
        filepath = "resources/{}".format(filename)
        urllib.request.urlretrieve(bed_file, filepath)
        return filepath
    else:
        return bed_file

# Check bed file for URL    
bed_file=check_bed_for_URL(snakemake.input[0])

# Function to create a bedpe file from a bed file
bed_list = []
with open(bed_file) as f:
    line = f.readlines()
    for name in line:
        bed_list.append(name.split())
df_bed = pd.DataFrame(
    bed_list, columns=["chrom", "start", "end", "name", "score", "strand"]
)
df_sense = df_bed.loc[df_bed["strand"] == "+"]
df_antisense = df_bed.loc[df_bed["strand"] == "-"]
# The dataframes for the sense and antisense strands need to be set to the same index so they can be merged again for the bedpe file
df_sense.reset_index(inplace=True)
df_antisense.reset_index(inplace=True)
data = [
    df_sense["chrom"],
    df_sense["start"],
    df_sense["end"],
    df_antisense["chrom"],
    df_antisense["start"],
    df_antisense["end"],
]
headers = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
df_bedpe = pd.concat(data, axis=1, keys=headers)
df_bedpe.to_csv(snakemake.output[0], header=None, index=None, sep="\t", mode="a")
