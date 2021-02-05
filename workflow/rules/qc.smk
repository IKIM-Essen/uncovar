rule fastqc:
    input:
        get_fastqs,
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip",
    log:
        "logs/fastqc/{sample}.log",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"


rule samtools_flagstat:
    input:
        "results/recal/ref~main/{sample}.bam",
    output:
        "results/qc/samtools_flagstat/{sample}.bam.flagstat",
    log:
        "logs/samtools/{sample}_flagstat.log",
    wrapper:
        "0.70.0/bio/samtools/flagstat"


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_fastqc.zip", sample=get_samples()),
        expand(
            "results/qc/samtools_flagstat/{sample}.bam.flagstat", sample=get_samples()
        ),
    output:
        "results/qc/multiqc.html",
    params:
        "",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "0.69.0/bio/multiqc"


# TODO Alexander and Thomas: add rules to detect contamination and perform QC
