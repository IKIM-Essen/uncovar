rule vcf_to_fasta:
    input:
        bcf="results/{date}/calls/ref~main/{sample}.bcf",
        fasta="resources/genomes/main.fasta",
        fai="resources/genomes/main.fasta.fai",
    output:
        "results/{date}/pseudoassembled-contigs/{sample}.fasta",
    params:
        min_prob_apply=config["assembly"]["min-variant-prob"],
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/vcf-to-fasta.py"
