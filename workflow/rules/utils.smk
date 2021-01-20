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
        "0.69.0/bio/tabix"


rule bam_index:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    log:
        "logs/bam-index/{prefix}.log",
    wrapper:
        "0.59.2/bio/samtools/index"
