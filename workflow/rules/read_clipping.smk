rule sort_bam:
    input:
        expand(
            "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam",
            ref=config["adapters"]["amplicon-reference"],
        ),
    output:
        temp("results/{date}/clipped-reads/{sample}.bam"),
    params:
        extra="-m 4G",
        tmp_dir="/tmp/",
    log:
        "logs/{date}/sort-bam/{sample}.log",
    threads: 8
    wrapper:
        "0.74.0/bio/samtools/sort"


rule clip_primer:
    input:
        sortbam="results/{date}/clipped-reads/{sample}.bam",
        sortindex="results/{date}/clipped-reads/{sample}.bam.bai",
        bed=config["adapters"]["amplicon-primers"],
        ref_fasta="resources/genomes/{reference}.fasta".format(
            reference=config["adapters"]["amplicon-reference"]
        ),
    output:
        clippedbam=temp("results/{date}/clipped-reads/{sample}.primerclipped.bam"),
        hardclippedbam=temp(
            "results/{date}/clipped-reads/{sample}.primerclipped.hard.bam"
        ),
        sorthardclippedbam=temp(
            "results/{date}/clipped-reads/{sample}.primerclipped.hard.namesorted.bam"
        ),
        fq1="results/{date}/clipped-reads/{sample}.1.fastq.gz",
        fq2="results/{date}/clipped-reads/{sample}.2.fastq.gz",
    params:
        dir=lambda w, input: os.path.dirname(input.sortbam),
        bam=lambda w, input: input.sortbam.split("/")[-1],
        dir_depth=lambda w, input: "".join(
            ["../"] * (len(input.sortbam.split("/")) - 1)
        ),
    log:
        "logs/{date}/primer-clipping/{sample}.log",
    conda:
        "../envs/bamclipper.yaml"
    threads: 4
    shell:
        """
        cd {params.dir}
        bamclipper.sh -b {params.bam} -p {params.dir_depth}{input.bed} -n {threads} > {params.dir_depth}{log} 2>&1
        cd {params.dir_depth}
        fgbio --sam-validation-stringency=LENIENT ClipBam -i {output.clippedbam} -o {output.hardclippedbam} -H true -r {input.ref_fasta} >> {log} 2>&1
        samtools sort  -@ {threads} -n {output.hardclippedbam} -o {output.sorthardclippedbam}  >> {log} 2>&1
        samtools fastq -@ {threads} {output.sorthardclippedbam} -1 {output.fq1} -2 {output.fq2}  >> {log} 2>&1
        """


rule sort_aln_for_plots:
    input:
        "results/{date}/clipped-reads/{sample}.primerclipped.hard.bam",
    output:
        "results/{date}/clipped-reads/{sample}.primerclipped.hard.sorted.bam",
    params:
        extra="-m 4G",
    log:
        "logs/{date}/sort-aln-for-plots/{sample}.log",
    wrapper:
        "0.77.0/bio/samtools/sort"


rule plot_primer_clipping:
    input:
        unclipped=expand_samples_for_date_amplicon(
            "results/{{date}}/clipped-reads/{sample}.bam"
        ),
        index_unclipped=expand_samples_for_date_amplicon(
            "results/{{date}}/clipped-reads/{sample}.bam.bai"
        ),
        clipped=expand_samples_for_date_amplicon(
            "results/{{date}}/clipped-reads/{sample}.primerclipped.hard.sorted.bam"
        ),
        index_clipped=expand_samples_for_date_amplicon(
            "results/{{date}}/clipped-reads/{sample}.primerclipped.hard.sorted.bam.bai"
        ),
    output:
        plot=report(
            "results/{date}/plots/primer-clipping-intervals.svg",
            caption="../report/amplicon-primer-clipping.rst",
            category="3. Sequencing Details",
            subcategory="4. Check for correct amplicon primer clipping",
        ),
    params:
        samples=lambda wildcards: get_samples_for_date(wildcards.date),
        bed=config["adapters"]["amplicon-primers"],
    log:
        "logs/{date}/plot-primer-clipping.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-primer-clipping.py"
