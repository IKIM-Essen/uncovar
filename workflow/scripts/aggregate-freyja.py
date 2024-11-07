# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import pandas as pd
import numpy as np
from natsort import index_natsorted

data = pd.DataFrame()
data2 = pd.DataFrame()

for file, sample, location, timestamp in zip(snakemake.input.demix, snakemake.params.sample, snakemake.params.location, snakemake.params.timestamp):
    lineages = []
    abundances = []
    # name, num = sample.split("-")
    with open(file, "r") as infile:
        for line in infile.read().splitlines():
            if line.startswith("lineages"):
                lineages = line.split("\t")[1].split(" ")
            elif line.startswith("abundances"):
                abundances = line.split("\t")[1].split(" ")
                rounded_abundances = [ '%.2f' % float(elem) for elem in abundances ]
        data.at[timestamp, location] = ", ".join([lin + ":" + ab for lin, ab in zip(lineages, rounded_abundances)])
        data2.at[timestamp, location] = len(lineages)

data.sort_index(key= lambda x: np.argsort(index_natsorted(data.index)), inplace=True)
data.to_csv(snakemake.output.all)
data2.to_csv(snakemake.output.all_count)


data = pd.DataFrame()

for file, sample in zip(snakemake.input.demix, snakemake.params.sample):
    lineages = []
    abundances = []
    with open(file, "r") as infile:
        for line in infile.read().splitlines():
            if line.startswith("lineages"):
                lineages = line.split("\t")[1].split(" ")
            elif line.startswith("abundances"):
                abundances = line.split("\t")[1].split(" ")
                rounded_abundances = [ '%.2f' % float(elem) for elem in abundances ]
        for lin, ab in zip(lineages, rounded_abundances):
            data.at[lin, sample] = ab

data.sort_index(key= lambda x: np.argsort(index_natsorted(data.index)), inplace=True)
# data.sort_index(key= lambda x: np.argsort(data.columns.to_series().str[1:].astype(int)), inplace=True)
data = data[sorted(data.columns, key=lambda x: tuple(map(int,x[-2:])))]
data.to_csv(snakemake.output.pivot)