rule fastqc:
    input:
        get_fastqs,
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip",
    log:
        "logs/fastqc/{sample}.log",
    threads: 1
    wrapper:
        "0.69.0/bio/fastqc"


# TODO include kraken after filtering, Kallisto
rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_fastqc.zip", sample=get_samples()),
        expand(
            "results/species-diversity/{sample}/{sample}.kreport2", sample=get_samples()
        ),
        expand("results/trimmed/{sample}.fastp.json", sample=get_samples()),
        expand("results/quast/{sample}/report.tsv", sample=get_samples()),
    output:
        "results/qc/multiqc.html",
    params:
        "--config config/multiqc_config.yaml",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "0.69.0/bio/multiqc"


# Analysis of species diversity present
rule species_diversity_before:
    input:
        db="resources/minikraken-8GB",
        reads=expand("results/trimmed/{{sample}}.{read}.fastq.gz", read=[1, 2]),
    output:
        classified_reads=expand(
            "results/species-diversity/{{sample}}/{{sample}}_{read}.classified.fasta",
            read=[1, 2],
        ),
        unclassified_reads=expand(
            "results/species-diversity/{{sample}}/{{sample}}_{read}.unclassified.fasta",
            read=[1, 2],
        ),
        kraken_output="results/species-diversity/{sample}/{sample}.kraken",
        report="results/species-diversity/{sample}/{sample}.kreport2",
    log:
        "logs/kraken/{sample}.log",
    params:
        classified=lambda w, output: "#".join(
            output.classified_reads[0].rsplit("_1", 1)
        ),
        unclassified=lambda w, output: "#".join(
            output.unclassified_reads[0].rsplit("_1", 1)
        ),
    threads: 8
    conda:
        "../envs/kraken.yaml"
    shell:
        "(kraken2 --db {input.db} --threads {threads} --unclassified-out {params.unclassified} --classified-out {params.classified} --report {output.report} --gzip-compressed --paired {input.reads} > {output.kraken_output}) 2> {log}"


# Plot Korna charts
rule create_krona_chart:
    input:
        kraken_output="results/species-diversity/{sample}/{sample}.kreport2",
        taxonomy_database="resources/krona/",
    output:
        "results/species-diversity/{sample}/{sample}.html",
    log:
        "logs/krona/{sample}.log",
    conda:
        "../envs/kraken.yaml"
    shell:
        "ktImportTaxonomy -m 3 -t 5 -tax {input.taxonomy_database} -o {output} {input} 2> {log}"


rule align_against_human:
    input:
        "resources/genomes/human-genome.fna.gz",
        expand("results/trimmed/{{sample}}.{read}.fastq.gz", read=[1, 2]),
    output:
        "results/ordered-contigs-human/{sample}.bam",
    log:
        "logs/minimap2/{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax asm5 {input} -o {output} 2> {log}"


rule extract_unmapped:
    input:
        "results/ordered-contigs-human/{sample}.bam",
    output:
        filtered="results/ordered-contigs-nonhuman/{sample}.bam",
        t=temp("results/ordered-contigs-nonhuman/{sample}/temp1.bam"),
        t2=temp("results/ordered-contigs-nonhuman/{sample}/temp2.bam"),
        t3=temp("results/ordered-contigs-nonhuman/{sample}/temp3.bam"),
    log:
        "logs/filter_human/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -u -f 4 -F 264 {input} > {output.t1}
        samtools view -u -f 8 -F 260 {input} > {output.t2}
        samtools view -u -f 12 -F 256 {input} > {output.t3}
        samtools merge -u - {output.t1} {output.t2} {output.t3} | samtools sort -n -o {output.filtered}
        """


rule merge_unmapped:
    


rule bamToFastq:
    input:
        "results/ordered-contigs-nonhuman/{sample}.bam",
    output:
        fq1="results/ordered-contigs-nonhuman/{sample}.1.fastq",
        fq2="results/ordered-contigs-nonhuman/{sample}.2.fastq",
    log:
        "logs/bamToFastq/{sample}.log",
    conda:
        "../envs/bedtools.yaml"
    shell:
        "bamtofastq -i {input} -fq {output.fq1} -fq2 {output.fq2}"


rule species_diversity_after:
    input:
        db="resources/minikraken-8GB",
        reads=expand(
            "results/ordered-contigs-nonhuman/{{sample}}.{read}.fastq", read=[1, 2]
        ),
    output:
        kraken_output="results/species-diversity-nonhuman/{sample}/{sample}.kraken",
        report="results/species-diversity-nonhuman/{sample}/{sample}.kreport2",
    log:
        "logs/kraken/{sample}_nonhuman.log",
    params:
    threads: 8
    conda:
        "../envs/kraken.yaml"
    shell:
        "(kraken2 --db {input.db} --threads {threads} --report {output.report} --gzip-compressed --paired {input.reads} > {output.kraken_output}) 2> {log}"


rule create_krona_chart_after:
    input:
        kraken_output="results/species-diversity-nonhuman/{sample}/{sample}.kreport2",
        taxonomy_database="resources/krona/",
    output:
        "results/species-diversity-nonhuman/{sample}/{sample}.html",
    log:
        "logs/krona/{sample}_nonhuman.log",
    conda:
        "../envs/kraken.yaml"
    shell:
        "ktImportTaxonomy -m 3 -t 5 -tax {input.taxonomy_database} -o {output} {input} 2> {log}"


# TODO Alexander and Thomas: add rules to detect contamination and perform QC
