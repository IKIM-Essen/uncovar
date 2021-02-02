rule tabix_gff:
    input:
        "{prefix}.{fmt}.gz",
    output:
        "{prefix}.{fmt}.gz.tbi",
    params:
        "-p {fmt}",
    log:
        "logs/tabix-{fmt}/{prefix}.log",
    wrapper:
        "0.70.0/bio/tabix"


rule bam_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/bam-index/{prefix}.log",
    wrapper:
        "0.70.0/bio/samtools/index"


rule faidx:
    input:
        "{prefix}.fasta",
    output:
        "{prefix}.fasta.fai",
    log:
        "logs/faidx/{prefix}.log",
    wrapper:
        "0.70.0/bio/samtools/faidx"