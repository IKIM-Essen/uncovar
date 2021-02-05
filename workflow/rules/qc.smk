# visualize quality of reads before trimming
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

# aggregate different qc outputs for all samples used
rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_fastqc.zip", sample=get_samples()),
    output:
        "results/qc/multiqc.html",
    params:
        "",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "0.69.0/bio/multiqc"


# TODO Alexander and Thomas: add rules to detect contamination and perform QC
