__author__ = "Johannes Köster"
__copyright__ = "Copyright 2021, Johannes Köster"
__email__ = "johannes.koester@uni-due.de"
__license__ = "MIT"

import os

from snakemake.shell import shell

pipeline = snakemake.params.pipeline
revision = snakemake.params.get("revision")
qs = snakemake.params.get("qs", [])
profile = snakemake.params.get("profile", [])
flags = snakemake.params.get("flags", [])

args = []


if isinstance(profile, str):
    profile = [profile]

if qs:
    args += ["-qs", str(qs)]

if revision:
    args += ["-revision", revision]

if profile:
    args += ["-profile", ",".join(profile)]

add_parameter = lambda name, value: args.append("--{} {}".format(name, value))

for name, files in snakemake.input.items():
    if isinstance(files, list):
        files = ",".join(files)
    add_parameter(name, files)
for name, value in snakemake.params.items():
    if (
        name != "pipeline"
        and name != "revision"
        and name != "profile"
        and name != "flags"
        and name != "qs"
    ):
        add_parameter(name, value)

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

args = " ".join(args)

if flags:
    args = " ".join([args, flags])

shell("nextflow run {pipeline} {args} {log}")
