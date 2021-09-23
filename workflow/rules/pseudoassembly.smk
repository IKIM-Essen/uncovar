rule vcf_to_fasta:
    input:
        bcf="results/{date}/calls/ref~main/{sample}.bcf",
        bam="results/{date}/recal/ref~main/{sample}.bam",
        bai="results/{date}/recal/ref~main/{sample}.bam.bai",
        fasta="resources/genomes/main.fasta",
        fai="resources/genomes/main.fasta.fai",
    output:
        "results/{date}/contigs/pseudoassembled/{sample}.fasta",
    params:
        min_prob_apply=config["assembly"]["min-variant-prob"],
        min_coverage=get_min_coverage,
    log:
        "logs/{date}/vcf-to-fasta/{sample}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/vcf-to-fasta.py"


rule compare_assemblies:
    input:
        assembly="results/{date}/contigs/polished/{sample}.fasta",
        pseudoassembly="results/{date}/contigs/pseudoassembled/{sample}.fasta",
    output:
        "results/{date}/aligned/assemblies/{sample}.bam",
    log:
        "logs/{date}/aligned/assemblies/{sample}log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 --eqx -ax asm5 {input.assembly} {input.pseudoassembly} -o {output} 2> {log}"


rule aggregate_assembly_comparisons:
    input:
        expand_samples_for_date("results/{{date}}/aligned/assemblies/{sample}.bam"),
    output:
        "results/{date}/tables/assembly_comparison.tsv",
    params:
        samples=lambda wildcards: get_samples_for_date(wildcards.date),
    log:
        "logs/{date}/aggregate-assembly-comparisons.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/aggregate-assembly-comparisons.py"
