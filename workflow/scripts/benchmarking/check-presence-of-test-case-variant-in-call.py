# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys

from pysam import VariantFile

# sys.stderr = open(snakemake.log[0], "w")


tests_cases = []

for path, pos, variant, test_case in zip(
    snakemake.input,
    snakemake.params.poses,
    snakemake.params.variants,
    snakemake.params.test_cases,
):
    with VariantFile(path, "rb") as variant_file:
        for record in variant_file.fetch():
            if int(pos) == record.pos:
                tests_cases.append(test_case)

with open(snakemake.output[0], "w") as f:
    for path in tests_cases:
        print(f"{path}", file=f)
