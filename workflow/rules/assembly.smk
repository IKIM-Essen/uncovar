rule assembly:
    input:
        fastq1="results/{date}/nonhuman-reads/{sample}.1.fastq.gz",
        fastq2="results/{date}/nonhuman-reads/{sample}.2.fastq.gz",
    output:
        temp("results/{date}/assembly/{sample}/{sample}.contigs.fa"),
    log:
        "logs/{date}/megahit/{sample}.log",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        sample=lambda wildcards: wildcards.sample,
    threads: 8
    conda:
        "../envs/megahit.yaml"
    shell:
        "(megahit -1 {input.fastq1} -2 {input.fastq2} --out-dir {params.outdir} -f && mv {params.outdir}/final.contigs.fa {output} ) 2> {log}"


rule align_unpolished_contigs:
    input:
        "resources/genomes/main.fasta",
        "results/{date}/assembly/{sample}/{sample}.contigs.fa",
    output:
        temp("results/{date}/aligned-unpolished-contigs/{sample}.bam"),
    log:
        "logs/{date}/minimap2/{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax asm5 {input} -o {output} 2> {log}"


rule quast_unpolished_contigs:
    input:
        reference="resources/genomes/main.fasta",
        bam="results/{date}/aligned-unpolished-contigs/{sample}.bam",
        fastas="results/{date}/assembly/{sample}/{sample}.contigs.fa",
    output:
        "results/{date}/quast-unpolished/{sample}/report.tsv",
    params:
        outdir=lambda x, output: os.path.dirname(output[0]),
    log:
        "logs/{date}/quast/{sample}.log",
    conda:
        "../envs/quast.yaml"
    threads: 8
    shell:
        "quast.py --min-contig 1 --threads {threads} -o {params.outdir} -r {input.reference} --bam {input.bam} {input.fastas} > {log} 2>&1"


rule order_contigs:
    input:
        contigs="results/{date}/assembly/{sample}/{sample}.contigs.fa",
        reference="resources/genomes/main.fasta",
    output:
        temp("results/{date}/ordered-contigs-all/{sample}.fasta"),
    log:
        "logs/{date}/ragoo/{sample}.log",
    params:
        outdir=lambda x, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/ragoo.yaml"
    shell:  # currently there is no conda package for mac available. Manuell download via https://github.com/malonge/RaGOO
        "(mkdir -p {params.outdir}/{wildcards.sample} && cd {params.outdir}/{wildcards.sample} && "
        "ragoo.py ../../../../{input.contigs} ../../../../{input.reference} && "
        "cd ../../../../ && mv {params.outdir}/{wildcards.sample}/ragoo_output/ragoo.fasta {output}) > {log} 2>&1"


rule filter_chr0:
    input:
        "results/{date}/ordered-contigs-all/{sample}.fasta",
    output:
        temp("results/{date}/ordered-contigs/{sample}.fasta"),
    log:
        "logs/{date}/ragoo/{sample}_cleaned.log",
    params:
        sample=lambda wildcards: wildcards.sample,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/ragoo_remove_chr0.py"


rule polish_contigs:
    input:
        fasta="results/{date}/ordered-contigs/{sample}.fasta",
        bcf="results/{date}/filtered-calls/ref~{sample}/{sample}.clonal.nofilter.bcf",
        bcfidx="results/{date}/filtered-calls/ref~{sample}/{sample}.clonal.nofilter.bcf.csi",
    output:
        report(
            "results/{date}/polished-contigs/{sample}.fasta",
            category="Assembly",
            caption="../report/assembly.rst",
        ),
    log:
        "logs/{date}/bcftools-consensus/{sample}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools consensus -f {input.fasta} {input.bcf} > {output} 2> {log}"


rule align_polished_contigs:
    input:
        "resources/genomes/main.fasta",
        "results/{date}/polished-contigs/{sample}.fasta",
    output:
        "results/{date}/aligned-polished-contigs/{sample}.bam",
    log:
        "logs/{date}/minimap2/{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax asm5 {input} -o {output} 2> {log}"


rule quast_polished_contigs:
    input:
        reference="resources/genomes/main.fasta",
        bam="results/{date}/aligned-polished-contigs/{sample}.bam",
        fastas="results/{date}/polished-contigs/{sample}.fasta",
    output:
        "results/{date}/quast-polished/{sample}/report.tsv",
    params:
        outdir=lambda x, output: os.path.dirname(output[0]),
    log:
        "logs/{date}/quast/{sample}.log",
    conda:
        "../envs/quast.yaml"
    threads: 8
    shell:
        "quast.py --min-contig 1 --threads {threads} -o {params.outdir} -r {input.reference} --bam {input.bam} {input.fastas} > {log} 2>&1"


# TODO blast smaller contigs to determine contamination?
