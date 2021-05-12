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


rule compare_assemblies:
    input:
        assembly="results/{date}/polished-contigs/{sample}.fasta",
        pseudoassembly="results/{date}/pseudoassembled-contigs/{sample}.fasta",
    output:
        "results/{date}/aligned/assemblies/{sample}.bam",
    log:
        "logs/{date}/aligned/assemblies/{sample}log"
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 --eqx -ax asm5 {input.assembly} {input.pseudoassembly} -o {output} 2> {log}"


rule aggregate_assembly_comparisons:
    input:
        lambda wildcards: expand(
            "results/{{date}}/aligned/assemblies/{sample}.bam",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        "results/{date}/tables/assembly_comparison.tsv",
    params:
        samples = lambda wildcards: get_samples_for_date(wildcards.date)
    log:
        "logs/{date}/aggregate-assembly-comparisons.log"
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/aggregate-assembly-comparisons.py"
