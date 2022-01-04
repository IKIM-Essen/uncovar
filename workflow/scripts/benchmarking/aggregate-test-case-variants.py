# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

all_variants = []
for path in snakemake.input:
    all_variants.append(pd.read_csv(path, sep="\t"))

pd.concat(all_variants).to_csv(snakemake.output[0], index=False, sep="\t")
