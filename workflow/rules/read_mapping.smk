rule bwa_index:
    input:
        get_reference(),
    output:
        multiext(
            "results/{date}/bwa/index/ref~{reference}.fasta",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    params:
        prefix=lambda w, output: os.path.splitext(output[0])[0],
    log:
        "logs/{date}/bwa-index/ref~{reference}.log",
    resources:
        mem_mb=369000,
    wrapper:
        "0.69.0/bio/bwa/index"


rule bwa_large_index:
    input:
        get_reference(),
    output:
        multiext(
            "resources/bwa/index/ref~{reference}.fasta",
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


rule map_reads:
    input:
        reads=get_reads,
        idx=get_bwa_index,
    output:
        temp("results/{date}/mapped/ref~{reference}/{sample}.bam"),
    log:
        "logs/{date}/bwa-mem/ref~{reference}/{sample}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra="",
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "0.69.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/{date}/mapped/ref~{reference}/{sample}.bam",
    output:
        bam=temp("results/{date}/dedup/ref~{reference}/{sample}.bam"),
        metrics="results/{date}/qc/dedup/ref~{reference}/{sample}.metrics.txt",
    log:
        "logs/{date}/picard/dedup/ref~{reference}/{sample}.log",
    params:
        "",
    wrapper:
        "0.69.0/bio/picard/markduplicates"


rule samtools_calmd:
    input:
        aln=get_recal_input,
        ref=get_reference(),
    output:
        "results/{date}/recal/ref~{reference}/{sample}.bam",
    log:
        "logs/{date}/samtools-calmd/ref~{reference}/{sample}.log",
    params:
        "-A",
    threads: 8
    wrapper:
        "0.69.0/bio/samtools/calmd"
