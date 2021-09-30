'# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.'
'# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)'
'# This file may not be copied, modified, or distributed'
'# except according to those terms.'
import sys
import pandas as pd
from snakemake.shell import shell

sys.stderr = open(snakemake.log[0], "w")

pangolin_results = pd.read_csv(snakemake.input.strain_call)
strain = pangolin_results.loc[0]["lineage"]

shell(
    "bcftools view -Ov {snakemake.input.bcfs} | (echo track name={snakemake.wildcards.target} description={strain}-{snakemake.wildcards.filter}; cat -) > {snakemake.output}"
)
