rule con_clipping:
    input:
        fastq1="results/{date}/nonhuman-reads/{sample}.1.fastq.gz",
        fastq2="results/{date}/nonhuman-reads/{sample}.2.fastq.gz",
        primers_left="/local/data/repos/snakemake-workflow-sars-cov2/resources/amp_primer_left.fasta",
        primers_right="/local/data/repos/snakemake-workflow-sars-cov2/resources/amp_primer_right.fasta",
    output:
        fastq1="results/{date}/con-clipped-reads/{sample}.1.fastq.gz",
        fastq2="results/{date}/con-clipped-reads/{sample}.2.fastq.gz",
    log:
        "logs/{date}/con-clipping/{sample}.log"
    threads: 8
    conda:
        "../envs/cutadapt.yaml"
    shell:
        "cutadapt -g file:{input.primers_left} -G file:{input.primers_right} -o {output.fastq1} -p {output.fastq2} {input.fastq1} {input.fastq2}"