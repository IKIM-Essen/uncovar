rule clipPrimer:
    input:
        bam="results/{date}/mapped/ref~MN908947/{sample}.bam",
        bed=config["adapters"]["amplicon-primers"],
    output:
        sortbam=temp("results/{date}/clipped-reads/{sample}.tmp.bam"),
        sortindex=temp("results/{date}/clipped-reads/{sample}.tmp.bam.bai"),
        tempbam=temp("results/{date}/clipped-reads/{sample}.tmp.primerclipped.bam"),
        sorttempbam=temp(
            "results/{date}/clipped-reads/{sample}.sort.tmp.primerclipped.bam"
        ),
        fq1="results/{date}/clipped-reads/{sample}.1.fastq.gz",
        fq2="results/{date}/clipped-reads/{sample}.2.fastq.gz",
    log:
        "logs/{date}/primer-clipping/{sample}.log",
    params:
        dir=lambda w, output: os.path.dirname(output.sortbam),
        bam=lambda w, output: output.sortbam.split("/")[-1],
        dir_depth=lambda w, output: "".join(
            ["../"] * (len(output.sortbam.split("/")) - 1)
        ),
    conda:
        "../envs/bamclipper.yaml"
    threads: 10
    shell:
        """
        samtools sort -@ {threads} -o {output.sortbam} {input.bam} > {log} 2>&1
        samtools index {output.sortbam} >> {log} 2>&1
        cd {params.dir}
        bamclipper.sh -b {params.bam} -p {params.dir_depth}{input.bed} -n {threads} 
        cd {params.dir_depth}
        samtools sort  -@ {threads} -n {output.tempbam} -o {output.sorttempbam}  >> {log} 2>&1
        samtools fastq -@ {threads} {output.sorttempbam} -1 {output.fq1} -2 {output.fq2}  >> {log} 2>&1
        """
