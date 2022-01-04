# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule fastp_pe:
    input:
        sample=get_fastqs,
    output:
        trimmed=temp(
            [
                "results/{date}/trimmed/fastp-pe/{sample}.1.fastq.gz",
                "results/{date}/trimmed/fastp-pe/{sample}.2.fastq.gz",
            ]
        ),
        html="results/{date}/trimmed/fastp-pe/{sample}.html",
        json="results/{date}/trimmed/fastp-pe/{sample}.fastp.json",
    params:
        adapters=get_adapters,
        extra="--qualified_quality_phred {} ".format(
            config["quality-criteria"]["illumina"]["min-PHRED"]
        )
        + "--length_required {}".format(
            config["quality-criteria"]["illumina"]["min-length-reads"]
        ),
    log:
        "logs/{date}/fastp/fastp-pe/{sample}.log",
    threads: 2
    wrapper:
        "0.70.0/bio/fastp"


rule fastp_se:
    input:
        sample=get_fastqs,
    output:
        trimmed=temp("results/{date}/trimmed/fastp-se/{sample}.fastq.gz"),
        html=temp("results/{date}/trimmed/fastp-se/{sample}.html"),
        json=temp("results/{date}/trimmed/fastp-se/{sample}.fastp.json"),
    params:
        adapters=get_adapters,
        extra="--qualified_quality_phred {} ".format(
            config["quality-criteria"]["illumina"]["min-PHRED"]
        )
        + "--length_required {}".format(
            config["quality-criteria"]["illumina"]["min-length-reads"]
        ),
    log:
        "results/{date}/trimmed/fastp-se/{sample}.log",
    threads: 2
    wrapper:
        "0.80.2/bio/fastp"
