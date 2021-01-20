rule bwa_index:
    input:
        "resources/genome.fasta"
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log:
        "logs/bwa-index.log"
    resources:
        mem_mb=369000
    cache: True
    wrapper:
        "0.69.0/bio/bwa/index"


rule map_reads:
    input:
        reads=get_fastqs,
        idx=rules.bwa_index.output
    output:
        temp("results/mapped/{sample}.bam")
    log:
        "logs/bwa-mem/{sample}.log"
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra="",
        sort="samtools",
        sort_order="coordinate"
    threads: 8
    wrapper:
        "0.69.0/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/mapped/{sample}.bam"
    output:
        bam="results/dedup/{sample}.bam",
        metrics="results/qc/dedup/{sample}.metrics.txt"
    log:
        "logs/picard/dedup/{sample}.log"
    params:
        ""
    wrapper:
        "0.69.0/bio/picard/markduplicates"


rule samtools_calmd:
    input:
        aln="results/dedup/{sample}.bam",
        ref="resources/genome.fasta"
    output:
        "results/recal/{sample}.bam"
    params:
        "-A"
    threads: 8
    wrapper:
        "0.69.0/bio/samtools/calmd"


rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    log:
        "logs/bam-index/{prefix}.log"
    wrapper:
        "0.59.2/bio/samtools/index"