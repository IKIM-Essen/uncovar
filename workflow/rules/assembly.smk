rule megahit_assembly:
    input:
        fastq1="results/trimmed/{sample}.1.fastq.gz",
        fastq2="results/trimmed/{sample}.2.fastq.gz",
    output:
        directory("results/assembly/{sample}_assembly.out"),
    log:
        "logs/megahit/{sample}.log",
    threads: 8
    conda:
        "../envs/megahit.yaml"
    shell:
        "megahit -1 {input.fastq1} -2 {input.fastq2} -o {output}"


rule metaSPAdes_assembly:
    input:
        fastq1="results/trimmed/{sample}.1.fastq.gz",
        fastq2="results/trimmed/{sample}.2.fastq.gz",
    output:
        "results/assembly/metaSPAdes/{sample}/contigs.fasta",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    log:
        "logs/metaSPAdes/{sample}.log",
    threads: 4
    conda:
        "../envs/spades.yaml"
    shell:
        "(metaspades.py -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir} -t {threads}) 2> {log}"
