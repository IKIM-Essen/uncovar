# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

sys.stderr = open(snakemake.log[0], "w")

with open(snakemake.output[0], "w") as out:
    print(*snakemake.params.mixtures, sep="\n", file=out)
