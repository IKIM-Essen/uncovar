__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os

from snakemake.shell import shell

pipeline = snakemake.params.pipeline
revision = snakemake.params.get("revision")
profile = snakemake.params.get("profile", [])
flags = snakemake.params.get("flags", [])

args = []


if isinstance(profile, str):
    profile = [profile]

if revision:
    args += ["-revision", revision]

if profile:
    args += ["-profile", ",".join(profile)]

# TODO pass threads in case of single job
# TODO limit parallelism in case of pipeline
# TODO handle other resources

add_parameter = lambda name, value: args.append("--{} {}".format(name, value))

for name, files in snakemake.input.items():
    if isinstance(files, list):
        # TODO check how multiple input files under a single arg are usually passed to nextflow
        files = ",".join(files)
    add_parameter(name, files)
for name, value in snakemake.params.items():
    if (
        name != "pipeline"
        and name != "revision"
        and name != "profile"
        and name != "flags"
    ):
        add_parameter(name, value)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

args = " ".join(args)

if flags:
    args = " ".join([args, flags])

shell("nextflow run {pipeline} {args} {log}")
