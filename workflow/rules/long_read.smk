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
        "NanoFilt --length {params.min_length} --quality {params.min_PHRED} --maxlength 600 {input} > {output} 2> {log}"  # --maxlength 700


rule convert2fasta:
    input:
        "results/{date}/filtered/nanofilt/{sample}.fastq",
    output:
        "results/{date}/filtered/nanofilt/{sample}.fasta",
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk seq -A  {input} > {output}"


rule map_for_capping:
    input:
        reads="results/{date}/filtered/nanofilt/{sample}.fasta",
        reference="resources/genomes/main.fasta",
    output:
        "results/{date}/minimappings/coverage/{sample}.paf",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -x map-ont {input.reference} {input.reads} -o {output} --secondary=no"


rule cap_cov_amp:
    input:
        primer="/home/simon/uncovar/.tests/resources/nCoV-2019.primer.bed",
        mappings="results/{date}/minimappings/coverage/{sample}.paf",
        reads="results/{date}/filtered/nanofilt/{sample}.fastq",
    output:
        "results/{date}/normalize_reads/{sample}_cap.fasta",
    script:
        "../scripts/amp_covcap_sampler.py"


rule canu_correct:
    input:
        # "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta",
        "results/{date}/normalize_reads/{sample}_cap.fasta",
    output:
        "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
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
        oeaConcurrency={params.concurrency}
        )
        2> {log}
        """  # 2>&1


# Intermediate number of threads (4-8) achieve best speedup of a+btrimming.
# For large files 8 threads help accelerate some, small files are processed faster with 4 threads.
rule porechop_adapter_barcode_trimming:
    input:
        # "results/{date}/normalize_reads/{sample}_cap.fasta",
        "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
    output:
        temp("results/{date}/corrected/trimmed/{sample}.corr.trim.fasta"),
    #    "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta",
    conda:
        "../envs/porechop.yaml"
    log:
        "logs/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.log",
    threads: 8
    shell:
        "porechop -i {input} -o {output} -t {threads} -v 1 > {log} 2>&1"


# rule map_corr:
#     input:
#         reads="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
#         reference="resources/genomes/main.fasta",
#     output:
#         alignments="results/{date}/minimappings/corrected/{sample}.paf",
#         dcreads="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta",
#     conda:
#         "../envs/minimap2.yaml"
#     shell:
#         """
#         minimap2 -x map-ont {input.reference} {input.reads} -o {output.alignments} --secondary=no &&
#         gzip -d {input.reads}
#         """


# rule trim_primers_corrected:
#     input:
#         # reads=get_fastqs,
#         primer="/home/simon/uncovar/.tests/resources/nCoV-2019.primer.bed",
#         mappings="results/{date}/minimappings/corrected/{sample}.paf",
#         reads="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta",
#     output:
#         "results/{date}/corrected/{sample}/{sample}.correctedReads.primerclip.fasta",
#     script:
#         "../scripts/map_trim.py"


rule map_raw:
    input:
        reads="results/{date}/normalize_reads/{sample}_cap.fasta",
        reference="resources/genomes/main.fasta",
    output:
        alignments="results/{date}/minimappings/trim_raw/{sample}.paf",
        # dcreads="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta",
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -x map-ont {input.reference} {input.reads} -o {output.alignments} --secondary=no
        """


rule trim_primers_raw:
    input:
        # reads=get_fastqs,
        primer="/home/simon/uncovar/.tests/resources/nCoV-2019.primer.bed",
        mappings="results/{date}/minimappings/trim_raw/{sample}.paf",
        reads="results/{date}/normalize_reads/{sample}_cap.fasta",
    output:
        "results/{date}/raw/trimmed/{sample}.raw.clip.fasta",
    script:
        "../scripts/map_trim.py"


# rule medaka_consensus_reference:
use rule assembly_polishing_ont as medaka_consensus_reference with:
    input:
        # Don´t ever use corrected reads as input for medaka, it is supposed to polish with raw-reads!
        # fasta="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
        # fasta="results/{date}/filtered/nanofilt/{sample}.fasta",
        # fasta="results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta",
        fasta="results/{date}/raw/trimmed/{sample}.raw.clip.fasta",
        reference="resources/genomes/main.fasta",
    output:
        # "results/{date}/consensus/medaka/{sample}/{sample}.fasta",
        "results/{date}/consensus/medaka/{sample}/consensus.fasta",
        # "results/{date}/polishing/medaka/{sample}/{sample}.consensus.fasta",


# polish consensus
rule bcftools_consensus_ont:
    input:
        fasta="results/{date}/consensus/medaka/{sample}/consensus.fasta",
        # fasta="results/{date}/consensus/medaka/{sample}/{sample}.fasta",
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
