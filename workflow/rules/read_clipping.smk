rule clip_primer:
    input:
        bam=expand(
            "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam",
            ref=config["adapters"]["amplicon-reference"],
        ),
        bed=config["adapters"]["amplicon-primers"],
        ref_fasta="resources/genomes/{reference}.fasta".format(reference=config["adapters"]["amplicon-reference"]),
    output:
        sortbam=temp("results/{date}/clipped-reads/{sample}.bam"),
        sortindex=temp("results/{date}/clipped-reads/{sample}.bam.bai"),
        clippedbam=temp("results/{date}/clipped-reads/{sample}.primerclipped.bam"),
        sortclippedbam=temp(
            "results/{date}/clipped-reads/{sample}.sort.primerclipped.bam"
        ),
        hardclippedbam=temp(
            "results/{date}/clipped-reads/{sample}.sort.primerclipped.hard.bam"
        ),
        fq1_initial=temp("results/{date}/clipped-reads/{sample}_initial.1.fastq.gz"),
        fq2_initial=temp("results/{date}/clipped-reads/{sample}_initial.2.fastq.gz"),
        fq1_sorted=temp("results/{date}/clipped-reads/{sample}.1.fastq"),
        fq2_sorted=temp("results/{date}/clipped-reads/{sample}.2.fastq"),
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
        bamclipper.sh -b {params.bam} -p {params.dir_depth}{input.bed} -n {threads} >> {params.dir_depth}{log} 2>&1
        cd {params.dir_depth}
        fgbio --sam-validation-stringency=LENIENT ClipBam -i {output.clippedbam} -o {output.hardclippedbam} -H true -r {input.ref_fasta} >> {log} 2>&1
        samtools sort  -@ {threads} -n {output.hardclippedbam} -o {output.sortclippedbam}  >> {log} 2>&1
        
        samtools fastq -@ {threads} {output.hardclippedbam} -1 {output.fq1_initial} -2 {output.fq2_initial}  >> {log} 2>&1
        zcat {output.fq1_initial} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\n\" > {output.fq1_sorted}
        zcat {output.fq2_initial} | paste - - - - | sort -k1,1 -t \" \" | tr \"\t\" \"\n\" > {output.fq2_sorted}
        gzip -c {output.fq1_sorted} > {output.fq1}
        gzip -c {output.fq2_sorted} > {output.fq2}
        """