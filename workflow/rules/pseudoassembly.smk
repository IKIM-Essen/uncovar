rule vcf_to_fasta:
    input:
        bcf="results/{date}/calls/ref~main/{sample}.bcf",
        bam="results/{date}/recal/ref~main/{sample}.bam",
        bai="results/{date}/recal/ref~main/{sample}.bam.bai",
        fasta="resources/genomes/main.fasta",
        fai="resources/genomes/main.fasta.fai",
    output:
        "results/{date}/pseudoassembled-contigs/{sample}.fasta",
    params:
        min_prob_apply=config["assembly"]["min-variant-prob"],
        min_coverage=get_min_coverage,
        sample=lambda wildcards: wildcards.sample,
    log:
        "logs/{date}/vcf-to-fasta/{sample}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/vcf-to-fasta.py"
