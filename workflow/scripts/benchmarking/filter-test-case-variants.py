# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

found_variants = pd.read_csv(snakemake.input[0], sep="\t")

different_probs = found_variants.loc[
    (found_variants["prob_present_illumina"] != found_variants["prob_present_ont"])
    & (found_variants["prob_present_illumina"] > 0)
    & (found_variants["prob_present_ont"] > 0)
]

illumina_only = found_variants.loc[
    (found_variants["prob_present_illumina"].isna())
    & (found_variants["prob_present_ont"] > 0)
]

ont_only = found_variants.loc[
    (found_variants["prob_present_illumina"] > 0)
    & (found_variants["prob_present_ont"].isna() > 0)
]


different_probs.to_csv(snakemake.output.different_probs, sep="\t", index=False)
illumina_only.to_csv(snakemake.output.illumina_only, sep="\t", index=False)
ont_only.to_csv(snakemake.output.ont_only, sep="\t", index=False)
