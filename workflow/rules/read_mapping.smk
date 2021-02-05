# bwa index for alignments
rule bwa_index:
    input:
        get_reference(),
    output:
        multiext(
            "results/bwa/index/ref~{reference}.fasta",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    params:
        prefix=lambda w, output: os.path.splitext(output[0])[0],
    log:
        "logs/bwa-index/ref~{reference}.log",
    resources:
        mem_mb=369000,
    wrapper:
        "0.69.0/bio/bwa/index"


# bwa mem alignment (local; based on smith-waterman with random seeding)
rule map_reads:
    input:
        reads=expand("results/trimmed/{{sample}}.{read}.fastq.gz", read=[1, 2]),
        idx=rules.bwa_index.output,
    output:
        temp("results/mapped/ref~{reference}/{sample}.bam"),
    log:
        "logs/bwa-mem/ref~{reference}/{sample}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra="",
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "0.69.0/bio/bwa/mem"


# ?????????????
rule mark_duplicates:
    input:
        "results/mapped/ref~{reference}/{sample}.bam",
    output:
        bam="results/dedup/ref~{reference}/{sample}.bam",
        metrics="results/qc/dedup/ref~{reference}/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/ref~{reference}/{sample}.log",
    params:
        "",
    wrapper:
        "0.69.0/bio/picard/markduplicates"


# calculate md sum for alignment ????????????
rule samtools_calmd:
    input:
        aln="results/dedup/ref~{reference}/{sample}.bam",
        ref=get_reference(),
    output:
        "results/recal/ref~{reference}/{sample}.bam",
    log:
        "logs/samtools-calmd/ref~{reference}/{sample}.log",
    params:
        "-A",
    threads: 8
    wrapper:
        "0.69.0/bio/samtools/calmd"
