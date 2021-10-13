# rule nanoQC_pre_trim:
#     input:
#         "input_data/{sample}.fastq"
#     output:
#         directory("results/{sample}/qc/initial")
#     log:
#         "logs/{sample}/{sample}_initialQC.log"
#     conda:
#         "envs/nanoqc.yaml"
#     shell:
#         "nanoQC {input} -o {output}"


# rule nanoQC_post_trim:
#     input:
#         "results/{sample}/trim-filt/{sample}_tf.fastq"
#     output:
#         directory("results/{sample}/qc/post_trimfilt_qc")
#     log:
#         "logs/{sample}/{sample}_post_trimfilt_qc.log"
#     conda:
#         "envs/nanoqc.yaml"
#     shell:
#         "nanoQC {input} -o {output}"


rule porechop_adapter_barcode_trimming:
    input:
        get_fastqs,
    output:
        "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fastq",
    conda:
        "../envs/porechop.yaml"
    log:
        "logs/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.log",
    threads: 2
    shell:
        "porechop -i {input} -o {output} --threads {threads} -v 4 > {log} 2>&1"


rule customize_primer_porechop:
    input:
        "resources/ARTIC_v3_adapters.py",
    output:
        "results/tables/replacement_notice.txt",
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
        fastq_in="results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fastq",
        repl_flag="results/tables/replacement_notice.txt",
    output:
        "results/{date}/trimmed/porechop/primer_clipped/{sample}.fastq",
    conda:
        "../envs/primechop.yaml"
    log:
        "logs/{date}/trimmed/porechop/primer_clipped/{sample}.log",
    threads: 2
    shell:
        "porechop -i {input.fastq_in} -o {output} --no_split --end_size 35 --extra_end_trim 0 --threads {threads} -v 4 > {log} 2>&1"


rule nanofilt:
    input:
        "results/{date}/trimmed/porechop/primer_clipped/{sample}.fastq",
    output:
        "results/{date}/trimmed/nanofilt/{sample}.fastq",
    log:
        "logs/{date}/nanofilt/{sample}.log",
    params:
        min_length=config["quality-criteria"]["ont"]["min-length-reads"],
        min_PHRED=config["quality-criteria"]["ont"]["min-PHRED"],
    conda:
        "../envs/nanofilt.yaml"
    shell:
        "NanoFilt --length {params.min_length} --quality {params.min_PHRED} {input} > {output} 2> {log}"


rule canu_correct:
    input:
        "results/{date}/trimmed/nanofilt/{sample}.fastq",
    output:
        "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
    log:
        "logs/{date}/canu/assemble/{sample}.log",
    params:
        outdir=get_output_dir,
        min_length=config["quality-criteria"]["ont"]["min-length-reads"],
    conda:
        "../envs/canu.yaml"
    threads: 6
    shell:
        "(canu -correct -nanopore {input} -p {wildcards.sample} -d {params.outdir}"
        " genomeSize=30k corOverlapper=minimap utgOverlapper=minimap obtOverlapper=minimap"
        " minOverlapLength=10 minReadLength={params.min_length} corMMapMerSize=10 corOutCoverage=50000"
        " corMinCoverage=0 maxInputCoverage=20000) "
        " 2> {log}"


rule spades_assemble_se:
    input:
        "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz",
    output:
        "results/{date}/assembly/{sample}/spades_se/{sample}.contigs.fasta",
    log:
        "logs/{date}/spades/se/{sample}.log",
    conda:
        "../envs/spades.yaml"
    params:
        outdir=get_output_dir,
    threads: 8
    shell:
        "(spades.py --corona -s {input} -o {params.outdir} -t {threads} && "
        " mv {params.outdir}/scaffolds.fasta {output})"
        " > {log} 2>&1"


# polish reference
rule medaka:
    input:
        fastq="results/{date}/trimmed/nanofilt/{sample}.fastq",
        reference="resources/genomes/main.fasta",
    output:
        "results/{date}/polished/medaka/{sample}/consensus.fasta",
    log:
        "logs/{date}/medaka/{sample}.log",
    params:
        outdir=get_output_dir,
        model=config["assembly"]["medaka_model"],
    conda:
        "../envs/medaka.yaml"
    threads: 4
    shell:
        "medaka_consensus -i {input.fastq} -o {params.outdir} -d {input.reference} -t {threads} -m {params.model} 2> {log}"


# Questions & Notes:
# Adjusted Nanoflit -> rmv max length and added min lenght
# Added canu correct before assembly
# Need an better way to provide amplicon primers for clipping
# We only want primer clipped de novo assembly
