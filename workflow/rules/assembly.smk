rule assembly:
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
