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


rule nanofilt:
    input:
        # get_reads_by_stage,
        get_fastqs,
    output:
        # "results/{date}/trimmed/porechop/primer_clipped/{sample}.fastq",
        # temp("results/{date}/trimmed/nanofilt/{sample}.fastq"),
        "results/{date}/filtered/nanofilt/{sample}.fastq",
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
        # "results/{barcode}/trim-filt/{barcode}_tf.fastq"
        "results/{date}/filtered/nanofilt/{sample}.fastq",
    output:
        "results/{date}/filtered/nanofilt/{sample}.fasta",
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk seq -A  {input} > {output}"


rule map_for_capping:
    input:
        # reads="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
        reads="results/{date}/filtered/nanofilt/{sample}.fasta",
        reference="resources/genomes/main.fasta",
    output:
        "results/{date}/minimappings/{sample}.paf",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -x map-ont {input.reference} {input.reads} -o {output} --secondary=no"


rule cap_cov_amp:
    input:
        # reads=get_fastqs,
        primer="/home/simon/uncovar/.tests/resources/nCoV-2019.primer.bed",
        mappings="results/{date}/minimappings/{sample}.paf",
        reads="results/{date}/filtered/nanofilt/{sample}.fastq",
    output:
        "results/{date}/normalize_reads/{sample}_cap.fasta",
    script:
        "../scripts/amp_covcap_sampler.py"


# Intermediate number of threads (4-8) achieve best speedup of a+btrimming.
# For large files 8 threads help accelerate some, small files are processed faster with 4 threads.
rule porechop_adapter_barcode_trimming:
    input:
        "results/{date}/normalize_reads/{sample}_cap.fasta",
    output:
        # temp("results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta"),
        "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta",
    conda:
        "../envs/porechop.yaml"
    log:
        "logs/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.log",
    threads: 8
    shell:
        "porechop -i {input} -o {output} -t {threads} -v 1 > {log} 2>&1"


rule canu_correct:
    input:
        # "results/{date}/normalize_reads/{sample}_cap.fasta",
        # "results/{date}/filtered/nanofilt/{sample}.fasta",
        # "results/{date}/trimmed/porechop/primer_clipped/{sample}.fasta",
        "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta",
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
        """


rule map_for_trimming:
    input:
        # reads="results/{date}/filtered/nanofilt/{sample}.fasta",
        reads="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
        reference="resources/genomes/main.fasta",
    output:
        alignments="results/{date}/minimappings/trimming/{sample}_trim.paf",
        dcreads="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta",
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -x map-ont {input.reference} {input.reads} -o {output.alignments} --secondary=no &&
        gzip -d {input.reads}
        """


rule trim_primers_corrected:
    input:
        # reads=get_fastqs,
        primer="/home/simon/uncovar/.tests/resources/nCoV-2019.primer.bed",
        mappings="results/{date}/minimappings/trimming/{sample}_trim.paf",
        reads="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta",
    output:
        "results/{date}/corrected/{sample}/{sample}.correctedReads.primerclip.fasta",
    script:
        "../scripts/map_trim.py"


# rule customize_primer_porechop:
#     input:
#         get_artic_primer,
#     output:
#         "results/.indicators/replacement_notice.txt",
#     conda:
#         "../envs/primechop.yaml"
#     log:
#         "logs/customize_primer_porechop.log",
#     shell:
#         "(cp {input} $CONDA_PREFIX/lib/python3.6/site-packages/porechop/adapters.py && "
#         'echo "replaced adpaters in adapter.py file in $CONDA_PREFIX/lib/python3.6/site-packages/porechop/adapters.py with ARTICv3 primers" > {output}) '
#         "2> {log}"


# Using a low number of threads (2-4) speed up primer-trimming significantly (>2x), even for large files,
# presumably due to the much higher number of target-sequences for trimming as compared
# to barcode+adapter-trimming. However, using only one thread is again very slow.
# rule porechop_primer_trimming:
#     input:
#         # fastq_in="results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta",
#         fastq_in="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
#         repl_flag="results/.indicators/replacement_notice.txt",
#     output:
#         "results/{date}/trimmed/porechop/primer_clipped/{sample}.corr.abpclip.fasta",
#     conda:
#         "../envs/primechop.yaml"
#     log:
#         "logs/{date}/trimmed/porechop/primer_clipped/{sample}.log",
#     threads: 2
#     shell:
#         """
#         (porechop -i {input.fastq_in} -o {output} --no_split --end_size 35 --extra_end_trim 0 -t {threads} -v 1) 2> {log}
#         rm results/.indicators/replacement_notice.txt
#         """


# rule count_fastq_reads:
#     input:
#         get_reads_by_stage,
#     output:
#         temp("results/{date}/tables/fastq-read-counts/{stage}~{sample}.txt"),
#     log:
#         "logs/{date}/count_reads/{stage}~{sample}.log",
#     conda:
#         "../envs/unix.yaml"
#     shell:
#         "echo $(( $(cat {input} | wc -l ) / 4)) > {output} 2> {log}"


# # Intermediate number of threads (4-8) achieve best speedup of a+btrimming.
# # For large files 8 threads help accelerate some, small files are processed faster with 4 threads.
# rule porechop_adapter_barcode_trimming:
#     input:
#         "results/{date}/normalize_reads/{sample}_cap.fasta",
#     output:
#         # temp("results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta"),
#         "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta"
#     conda:
#         "../envs/porechop.yaml"
#     log:
#         "logs/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.log",
#     threads: 8
#     shell:
#         "porechop -i {input} -o {output} -t {threads} -v 1 > {log} 2>&1"


# rule customize_primer_porechop:
#     input:
#         get_artic_primer,
#     output:
#         "results/.indicators/replacement_notice.txt",
#     conda:
#         "../envs/primechop.yaml"
#     log:
#         "logs/customize_primer_porechop.log",
#     shell:
#         "(cp {input} $CONDA_PREFIX/lib/python3.6/site-packages/porechop/adapters.py && "
#         'echo "replaced adpaters in adapter.py file in $CONDA_PREFIX/lib/python3.6/site-packages/porechop/adapters.py with ARTICv3 primers" > {output}) '
#         "2> {log}"


# # Using a low number of threads (2-4) speed up primer-trimming significantly (>2x), even for large files,
# # presumably due to the much higher number of target-sequences for trimming as compared
# # to barcode+adapter-trimming. However, using only one thread is again very slow.
# rule porechop_primer_trimming:
#     input:
#         fastq_in=(
#             "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta"
#         ),
#         repl_flag="results/.indicators/replacement_notice.txt",
#     output:
#         "results/{date}/trimmed/porechop/primer_clipped/{sample}.fasta",
#     conda:
#         "../envs/primechop.yaml"
#     log:
#         "logs/{date}/trimmed/porechop/primer_clipped/{sample}.log",
#     threads: 2
#     shell:
#         """
#         (porechop -i {input.fastq_in} -o {output} --no_split --end_size 35 --extra_end_trim 0 -t {threads} -v 1) 2> {log}
#         rm results/.indicators/replacement_notice.txt
#         """


# rule medaka_consensus_reference:
use rule assembly_polishing_ont as medaka_consensus_reference with:
    input:
        # Don´t ever use corrected reads as input for medaka, it is supposed to polish with rwa-reads!
        # fasta="results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
        # fasta="results/{date}/filtered/nanofilt/{sample}.fasta",
        "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fasta",
        reference="resources/genomes/main.fasta",
    output:
        "results/{date}/consensus/medaka/{sample}/{sample}.fasta",


# polish consensus
rule bcftools_consensus_ont:
    input:
        fasta="results/{date}/consensus/medaka/{sample}/{sample}.fasta",
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
