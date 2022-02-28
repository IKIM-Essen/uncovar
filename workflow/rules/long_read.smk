# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule nanoQC:
    input:
        get_reads_by_stage,
    output:
        temp("results/{date}/qc/nanoQC/{sample}/{stage}/nanoQC.html"),
    log:
        "logs/{date}/nanoQC/{sample}_{stage}.log",
    params:
        outdir=get_output_dir,
    conda:
        "../envs/nanoqc.yaml"
    shell:
        "nanoQC {input} -o {params.outdir} > {log} 2>&1"


rule count_fastq_reads:
    input:
        get_reads_by_stage,
    output:
        temp("results/{date}/tables/fastq-read-counts/{stage}~{sample}.txt"),
    log:
        "logs/{date}/count_reads/{stage}~{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "echo $(( $(cat {input} | wc -l ) / 4)) > {output} 2> {log}"


rule nanofilt:
    input:
        get_fastqs,
    output:
        temp("results/{date}/filtered/nanofilt/{sample}.fastq"),
    log:
        "logs/{date}/nanofilt/{sample}.log",
    params:
        min_length=config["quality-criteria"]["ont"]["min-length-reads"],
        min_PHRED=config["quality-criteria"]["ont"]["min-PHRED"],
    conda:
        "../envs/nanofilt.yaml"
    shell:
        "NanoFilt --length {params.min_length} --quality {params.min_PHRED} --maxlength 700 {input} > {output} 2> {log}"


rule downsample_and_trim_raw:
    input:
        primer="/home/simon/uncovar/.tests/resources/nCoV-2019.primer.bed",
        reads="results/{date}/filtered/nanofilt/{sample}.fastq",
        ref_genome="resources/genomes/main.fasta",
    output:
        "results/{date}/norm_trim_raw_reads/{sample}/{sample}.cap.fasta",
        "results/{date}/norm_trim_raw_reads/{sample}/{sample}.cap.clip.fasta",
    params:
        outdir="results/{date}/norm_trim_raw_reads/{sample}",
    log:
        "results/{date}/norm_trim_raw_reads/{sample}/notramp.log",
    conda:
        "../envs/notramp.yaml"
    shell:
        "notramp -a -r {input.reads} -p {input.primer} -g {input.ref_genome} -o {params.outdir}"


# rule medaka_consensus_reference:
use rule assembly_polishing_ont as medaka_consensus_reference with:
    input:
        # Don´t ever use corrected reads as input for medaka, it is supposed to polish with raw-reads!
        fasta="results/{date}/norm_trim_raw_reads/{sample}/{sample}.cap.clip.fasta",
        reference="resources/genomes/main.fasta",
    output:
        "results/{date}/consensus/medaka/{sample}/consensus.fasta",


rule canu_correct:
    input:
        "results/{date}/norm_trim_raw_reads/{sample}/{sample}.cap.fasta",
    output:
        "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta",
    log:
        "logs/{date}/canu/correct/{sample}.log",
    params:
        outdir=get_output_dir,
        concurrency=lambda w, threads: get_canu_concurrency(threads),
        min_length=config["quality-criteria"]["ont"]["min-length-reads"],
        for_testing=lambda w, threads: get_if_testing(
            f"corThreads={threads} redMemory=6 oeaMemory=6"
        ),
    conda:
        "../envs/canu.yaml"
    threads: 16
    shell:
        """
        ( if [ -d {params.outdir} ]; then rm -Rf {params.outdir}; fi &&
        canu -correct -nanopore {input} -p {wildcards.sample} -d {params.outdir} genomeSize=30k minOverlapLength=10 minReadLength=200 \
        useGrid=false {params.for_testing} \
        corMMapMerSize=10 corOutCoverage=50000 corMinCoverage=0 maxInputCoverage=20000 \
        corOverlapper=minimap utgOverlapper=minimap obtOverlapper=minimap \
        corConcurrency={params.concurrency} \
        cormhapConcurrency={params.concurrency} cormhapThreads={params.concurrency} \
        cormmapConcurrency={params.concurrency} cormmapThreads={params.concurrency} \
        obtmmapConcurrency={params.concurrency} obtmmapThreads={params.concurrency} \
        utgmmapConcurrency={params.concurrency} utgmmapThreads={params.concurrency} \
        redConcurrency={params.concurrency} redThreads={params.concurrency} \
        ovbConcurrency={params.concurrency} \
        ovsConcurrency={params.concurrency} \
        oeaConcurrency={params.concurrency})
        gzip -d {output}.gz
        2> {log}
        """


rule clip_adbc_corrected:
    input:
        primer="/home/simon/uncovar/.tests/resources/nCoV-2019.primer.bed",
        reads="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta",
        ref_genome="resources/genomes/main.fasta",
    output:
        "results/{date}/norm_trim_corr_reads/{sample}/{sample}.clip.fasta",
    params:
        outdir="results/{date}/norm_trim_corr_reads/{sample}",
    log:
        "results/{date}/norm_trim_corr_reads/{sample}/notramp.log",
    conda:
        "../envs/notramp.yaml"
    shell:
        "notramp -t --incl_prim -r {input.reads} -p {input.primer} -g {input.ref_genome} -o {params.outdir}"


rule bcftools_consensus_ont:
    input:
        fasta="results/{date}/consensus/medaka/{sample}/consensus.fasta",
        bcf="results/{date}/filtered-calls/ref~{sample}/{sample}.subclonal.high+moderate-impact.bcf",  # clonal vs. subclonal?
        bcfidx="results/{date}/filtered-calls/ref~{sample}/{sample}.subclonal.high+moderate-impact.bcf.csi",
    output:
        "results/{date}/consensus/bcftools/{sample}.fasta",
    log:
        "logs/{date}/bcftools-consensus-ont/{sample}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools consensus -f {input.fasta} {input.bcf} > {output} 2> {log}"


rule rename_consensus:
    input:
        "results/{date}/consensus/bcftools/{sample}.fasta",
    output:
        report(
            "results/{date}/contigs/consensus/{sample}.fasta",
            category="4. Sequences",
            subcategory="3. Consensus Sequences",
            caption="../report/assembly_consensus.rst",
        ),
    log:
        "logs/{date}/rename-consensus-fasta/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(cp {input} {output} && "
        ' sed -i "1s/.*/>{wildcards.sample}/" {output})'
        " 2> {log}"
