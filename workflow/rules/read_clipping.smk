rule clip_primer:
    input:
        bam=expand(
            "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam",
            ref=config["adapters"]["amplicon-reference"],
        ),
        bed=config["adapters"]["amplicon-primers"],
        ref_fasta="resources/genomes/{reference}.fasta".format(
            reference=config["adapters"]["amplicon-reference"]
        ),
    output:
        sortbam=temp("results/{date}/clipped-reads/{sample}.bam"),
        clippedbam=temp("results/{date}/clipped-reads/{sample}.primerclipped.bam"),
        hardclippedbam=temp(
            "results/{date}/clipped-reads/{sample}.primerclipped.hard.bam"
        ),
        sorthardclippedbam=temp(
            "results/{date}/clipped-reads/{sample}.primerclipped.hard.sorted.bam"
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

        cd {params.dir}
        bamclipper.sh -b {params.bam} -p {params.dir_depth}{input.bed} -n {threads} >> {params.dir_depth}{log} 2>&1
        cd {params.dir_depth}
        fgbio --sam-validation-stringency=LENIENT ClipBam -i {output.clippedbam} -o {output.hardclippedbam} -H true -r {input.ref_fasta} >> {log} 2>&1
        samtools sort  -@ {threads} -n {output.hardclippedbam} -o {output.sorthardclippedbam}  >> {log} 2>&1
        samtools fastq -@ {threads} {output.sorthardclippedbam} -1 {output.fq1} -2 {output.fq2}  >> {log} 2>&1
        """


rule clp_aln:
    input:
        reads=[
            "results/{date}/clipped-reads/{sample}.1.fastq.gz",
            "results/{date}/clipped-reads/{sample}.2.fastq.gz",
        ],
        idx=get_bwa_index,
    output:
        "results/{date}/clipped-aln/ref~{reference}/{sample}.bam",
    log:
        "logs/{date}/clipped-aln/ref~{reference}/{sample}.log",
    params:
        index=lambda w, input: os.path.splitext(input.idx[0])[0],
        extra="",
        sort="samtools",
        sort_order="coordinate",
    threads: 8
    wrapper:
        "0.69.0/bio/bwa/mem"


rule sort_aln_for_plots:
    input:
        "results/{date}/clipped-reads/{sample}.primerclipped.hard.bam",
    output:
        "results/{date}/clipped-reads/{sample}.primerclipped.hard.c_sort.bam",
    log:
        "logs/{date}/sort-aln-for-plots/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -o {output} {input} > {log} 2>&1"


rule plot_primer_clipping:
    input:
        unclipped=lambda wildcards: expand(
            "results/{{date}}/clipped-reads/{sample}.bam",
            sample=get_samples_for_date(wildcards.date),
        ),
        index_unclipped=lambda wildcards: expand(
            "results/{{date}}/clipped-reads/{sample}.bam.bai",
            sample=get_samples_for_date(wildcards.date),
        ),
        clipped=lambda wildcards: expand(
            "results/{{date}}/clipped-reads/{sample}.primerclipped.hard.c_sort.bam",
            sample=get_samples_for_date(wildcards.date),
        ),
        index_clipped=lambda wildcards: expand(
            "results/{{date}}/clipped-reads/{sample}.primerclipped.hard.c_sort.bam.bai",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        plot="results/{date}/plots/all.svg",
    log:
        "logs/{date}/plot-primer-clipping/{date}.log",
    threads: 16
    params:
        samples=lambda wildcards: get_samples_for_date(wildcards.date),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-primer-clipping.py"
