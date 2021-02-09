rule assembly:
    input:
        fastq1="results/ordered-contigs-nonhuman/{sample}.1.fastq.gz",
        fastq2="results/ordered-contigs-nonhuman/{sample}.2.fastq.gz",
    output:
        "results/assembly/{sample}/{sample}.contigs.fa",
    log:
        "logs/megahit/{sample}.log",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        sample=lambda wildcards: wildcards.sample,
    threads: 8
    conda:
        "../envs/megahit.yaml"
    shell:
        "(megahit -1 {input.fastq1} -2 {input.fastq2} --out-dir {params.outdir} -f && mv {params.outdir}/{params.sample}.contigs.fa {output} ) 2> {log}"


rule align_contigs:
    input:
        "resources/genomes/main.fasta",
        "results/assembly/{sample}/{sample}.contigs.fa",
    output:
        "results/ordered-contigs/{sample}.bam",
    log:
        "logs/minimap2/{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax asm5 {input} -o {output} 2> {log}"


rule quast:
    input:
        reference="resources/genomes/main.fasta",
        bam="results/ordered-contigs/{sample}.bam",
        fastas="results/assembly/{sample}/{sample}.contigs.fa",
    output:
        "results/quast/{sample}/report.tsv",
    params:
        outdir=lambda x, output: os.path.dirname(output[0]),
    log:
        "logs/quast/{sample}.log",
    conda:
        "../envs/quast.yaml"
    threads: 8
    shell:
        "quast.py --threads {threads} -o {params.outdir} -r {input.reference} --bam {input.bam} {input.fastas} > {log} 2>&1"


rule order_contigs:
    input:
        contigs="results/assembly/{sample}/{sample}.contigs.fa",
        reference="resources/genomes/main.fasta",
    output:
        "results/ordered-contigs-all/{sample}.fasta",
    log:
        "logs/ragoo/{sample}.log",
    params:
        outdir=lambda x, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/ragoo.yaml"
    shell:  # currently there is no conda package for mac available. Manuell download via https://github.com/malonge/RaGOO
        "(mkdir -p {params.outdir}/{wildcards.sample} && cd {params.outdir}/{wildcards.sample} && "
        "ragoo.py ../../../{input.contigs} ../../../{input.reference} && "
        "cd ../../../ && mv {params.outdir}/{wildcards.sample}/ragoo_output/ragoo.fasta {output}) > {log} 2>&1"


rule filter_chr0:
    input:
        "results/ordered-contigs-all/{sample}.fasta",
    output:
        "results/ordered-contigs/{sample}.fasta",
    log:
        "logs/ragoo/{sample}_cleaned.log",
    threads: 8
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/ragoo_remove_chr0.py"


rule polish_contigs:
    input:
        fasta="results/ordered-contigs/{sample}.fasta",
        bcf="results/filtered-calls/ref~{sample}/{sample}.clonal.bcf",
        bcfidx="results/filtered-calls/ref~{sample}/{sample}.clonal.bcf.csi",
    output:
        report(
            "results/polished-contigs/{sample}.fasta",
            category="Assembly",
            caption="../report/assembly.rst",
        ),
    log:
        "logs/bcftools-consensus/{sample}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools consensus -f {input.fasta} {input.bcf} > {output} 2> {log}"


# TODO blast smaller contigs to determine contamination?
# TODO add Quast for polished contigs?
