import subprocess
import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

pangolin_results = pd.read_csv(snakemake.input.strain_call)
strain = pangolin_results.loc[0]["lineage"]

subprocess.call("bcftools view -Ov %s | (echo track name=%s description=%s-%s; cat -) > %s" % (snakemake.input.bcfs, snakemake.wildcards.target, strain, snakemake.wildcards.filter, snakemake.output), shell=True)