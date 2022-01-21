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

# Intermediate number of threads (4-8) achieve best speedup of a+btrimming.
# For large files 8 threads help accelerate some, small files are processed faster with 4 threads.  
rule porechop_adapter_barcode_trimming:
    input:
        get_fastqs,
    output:
        temp("results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fastq"),
    conda:
        "../envs/porechop.yaml"
    log:
        "logs/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.log",
    threads: 8
    shell:
        "porechop -i {input} -o {output} -t {threads} -v 1 > {log} 2>&1"

# Using a low number of threads (2-4) speed up primer-trimming significantly (>2x), even for large files,
# presumably due to the much higher number of target-sequences for trimming as compared 
# to barcode+adapter-trimming. However, using only one thread is again very slow.
rule customize_primer_porechop:
    input:
        get_artic_primer,
    output:
        "results/.indicators/replacement_notice.txt",
    conda:
        "../envs/primechop.yaml"
    log:
        "logs/customize_primer_porechop.log",
    shell:
        "(cp {input} $CONDA_PREFIX/lib/python3.6/site-packages/porechop/adapters.py && "
        'echo "replaced adpaters in adapter.py file in $CONDA_PREFIX/lib/python3.6/site-packages/porechop/adapters.py with ARTICv3 primers" > {output}) '
        "2> {log}"


rule porechop_primer_trimming:
    input:
        fastq_in=(
            "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fastq"
        ),
        repl_flag="results/.indicators/replacement_notice.txt",
    output:
        temp("results/{date}/trimmed/porechop/primer_clipped/{sample}.fastq"),
    conda:
        "../envs/primechop.yaml"
    log:
        "logs/{date}/trimmed/porechop/primer_clipped/{sample}.log",
    threads: 2
    shell:
        """
        (porechop -i {input.fastq_in} -o {output} --no_split --end_size 35 --extra_end_trim 0 -t {threads} -v 1) 2> {log}
        rm results/.indicators/replacement_notice.txt
        """


rule nanofilt:
    input:
        "results/{date}/trimmed/porechop/primer_clipped/{sample}.fastq",
    output:
        temp("results/{date}/trimmed/nanofilt/{sample}.fastq"),
    log:
        "logs/{date}/nanofilt/{sample}.log",
    params:
        min_length=config["quality-criteria"]["ont"]["min-length-reads"],
        min_PHRED=config["quality-criteria"]["ont"]["min-PHRED"],
    conda:
        "../envs/nanofilt.yaml"
    shell:
        "NanoFilt --length {params.min_length} --quality {params.min_PHRED} --maxlength 500 {input} > {output} 2> {log}"


rule canu_correct:
    input:
        "results/{date}/trimmed/nanofilt/{sample}.fastq",
    output:
        "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
    log:
        "logs/{date}/canu/assemble/{sample}.log",
    params:
        outdir=get_output_dir,
        concurrency = lambda w, threads: int(threads/ 4) ,
        min_length=config["quality-criteria"]["ont"]["min-length-reads"],

        for_testing=lambda w, threads: get_if_testing(
            f"corThreads={threads} redThreads={threads} redMemory=6 oeaMemory=6"
        ),
    conda:
        "../envs/canu.yaml"
    threads: 16
    shell:
        """
        ( if [ -d {params.outdir} ]; then rm -Rf {params.outdir}; fi &&
        canu -correct -nanopore {input} -p {wildcards.sample} -d {params.outdir} genomeSize=30k minOverlapLength=10 minReadLength=200 \
        useGrid=false \
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
        oeaConcurrency={params.concurrency}
        )
        2> {log}
        """


# rule medaka_consensus_reference:
use rule assembly_polishing_ont as medaka_consensus_reference with:
    input:
        fasta="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
        reference="resources/genomes/main.fasta",
    output:
        temp("results/{date}/consensus/medaka/{sample}/{sample}.fasta"),


# polish consensus
rule bcftools_consensus_ont:
    input:
        fasta="results/{date}/consensus/medaka/{sample}/{sample}.fasta",
        bcf="results/{date}/filtered-calls/ref~{sample}/{sample}.subclonal.high+moderate-impact.bcf",  # clonal vs. subclonal?
        bcfidx="results/{date}/filtered-calls/ref~{sample}/{sample}.subclonal.high+moderate-impact.bcf.csi",
    output:
        temp("results/{date}/consensus/bcftools/{sample}.fasta"),
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
