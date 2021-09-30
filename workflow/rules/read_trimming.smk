'# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.'
'# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)'
'# This file may not be copied, modified, or distributed'
'# except according to those terms.'
rule fastp_pe:
    input:
        sample=get_fastqs,
    output:
        trimmed=temp(
            [
                "results/{date}/trimmed/{sample}.1.fastq.gz",
                "results/{date}/trimmed/{sample}.2.fastq.gz",
            ]
        ),
        html="results/{date}/trimmed/{sample}.html",
        json="results/{date}/trimmed/{sample}.fastp.json",
    params:
        adapters=get_adapters,
        extra="--qualified_quality_phred {} ".format(
            config["RKI-quality-criteria"]["illumina"]["min-PHRED"]
        )
        + "--length_required {}".format(
            config["RKI-quality-criteria"]["illumina"]["min-length-reads"]
        ),
    log:
        "logs/{date}/fastp/{sample}.log",
    threads: 2
    wrapper:
        "0.70.0/bio/fastp"
