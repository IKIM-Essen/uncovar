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
        get_fastqs
    output:
        "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fastq"
    conda:
        "../envs/porechop.yaml"
    log:
        "logs/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.log",
    threads: 2
    shell:
        "porechop -i {input} -o {output} --threads {threads} -v 4 > {log} 2>&1"


# rule customize_primer_porechop:
#     output:
#         "workflow/report/replacement_notice.txt"
#     conda:
#         "envs/primechop.yaml"
#     shell:
#         #rm $CONDA_PREFIX/lib/python3.6/site-packages/porechop/adapters.py
#         """
#         cp ./resources/ARTIC_v3_adapters.py $CONDA_PREFIX/lib/python3.6/site-packages/porechop/adapters.py
#         echo "replaced adpaters in adapter.py file with ARTICv3 primers" > {output}
#         """


# rule porechop_primer_trimming:
#     input:
#         fastq_in="results/{sample}_abtrim/{sample}.fastq",
#         repl_flag="workflow/report/replacement_notice.txt"
#     output:
#         "results/{sample}_primertrim/{sample}.fastq"
#     conda:
#         "envs/primechop.yaml"
#     threads: 2
#     shell:
#         "porechop -i {input.fastq_in} -o {output} --no_split --end_size 35 --extra_end_trim 0 --threads {threads} -v 0"


rule nanofilt:
    input:
        "results/{date}/trimmed/porechop/adapter_barcode_trimming/{sample}.fastq"
    output:
        "results/{date}/trimmed/nanofilt/{sample}.fastq"
    log:
        "logs/{date}/nanofilt/{sample}.log"
    params:
        min_length=config["RKI-quality-criteria"]["ont"]["min-length-reads"],
        min_PHRED=config["RKI-quality-criteria"]["ont"]["min-PHRED"]
    conda:
        "../envs/nanofilt.yaml"
    shell:
        "NanoFilt --length {params.min_length} --quality {params.min_PHRED} {input} > {output} 2> {log}"


# rule seqtk_convert2fasta:
#     input:
#         "results/{date}/trimmed/nanofilt/{sample}.fastq"
#     output:
#         "results/{date}/trimmed/nanofilt_fasta/{sample}.fasta"
#     log:
#         "logs/{date}/nanofilt/{sample}.log"
#     conda:
#         "../envs/seqtk.yaml"
#     shell:
#         "seqtk seq -A {input} > {output} 2> {log}"


rule canu_correct:
    input:
        "results/{date}/trimmed/nanofilt/{sample}.fastq"
    output:
        "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz"
    log:
        "logs/{date}/canu/assemble/{sample}.log"
    params:
        outdir=get_output_dir
    conda:
        "../envs/canu.yaml"
    threads: 6
    shell:
        "(canu -correct -nanopore {input} -p {wildcards.sample} -d {params.outdir} "
        "genomeSize=30k corOverlapper=minimap utgOverlapper=minimap obtOverlapper=minimap "
        "minOverlapLength=10 minReadLength=200 corMMapMerSize=10 corOutCoverage=50000 "
        "corMinCoverage=0 maxInputCoverage=20000) "
        "2> {log}"


rule assembly_canu:
    input:
        "results/{date}/corrected/{sample}/{sample}.correctedReads.fasta.gz"
    output:
        "results/{date}/assembly/{sample}/canu/{sample}.contigs.fasta"
    params:
        outdir=get_output_dir
    log:
        "logs/{date}/canu/assemble/{sample}.log"
    conda:
        "../envs/canu.yaml"
    threads: 6
    shell:
        "(canu -assemble -nanopore -corrected {input} -p {wildcards.sample} -d {params.outdir} genomeSize=30k corOverlapper=minimap "
        "utgOverlapper=minimap obtOverlapper=minimap minOverlapLength=20 minReadLength=300 corMMapMerSize=10) "
        "2> {log}"


# polish reference
rule medaka:
    input:
        fastq="results/{date}/trimmed/nanofilt/{sample}.fastq",
        reference="resources/genomes/main.fasta",
    output:
        "results/{date}/polished/medaka/{sample}/consensus.fasta"
    log:
        "logs/{date}/medaka/{sample}.log"
    params:
        outdir=get_output_dir,
        model = config["assembly"]["medaka_model"]
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