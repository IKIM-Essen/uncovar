import tempfile
import shutil
from pathlib import Path

from snakemake.shell import shell

with open(snakemake.input.depth) as infile:
    depth = int(infile.read())

# https://bitbucket.org/jbaaijens/virus-vg/src/master/
min_node_abundance = depth * 0.005
min_strain_abundance = depth * 0.01

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

with tempfile.TemporaryDirectory() as tmp:
    shell(
        "cd {tmp}; "
        "{snakemake.input.virus_vg}/scripts/optimize_strains.py "
        "-m {min_node_abundance} -c {min_strain_abundance} "
        "{snakemake.input.abundances} {snakemake.input.gfa} 2> {log}"
    )
    tmp = Path(tmp)

    shutil.move(tmp / "haps.final.fasta", snakemake.output[0])