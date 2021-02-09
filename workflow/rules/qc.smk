rule fastqc_R1:
    input:
        lambda wildcards: pep.sample_table.loc[wildcards.sample][["fq1"]]
    output:
        html="results/qc/fastqc/{sample}.1.html",
        zip="results/qc/fastqc/{sample}.1_fastqc.zip",
    log:
        "logs/fastqc/{sample}.1.log",
    threads: 8
    wrapper:
        "0.69.0/bio/fastqc"


rule fastqc_R2:
    input:
        lambda wildcards: pep.sample_table.loc[wildcards.sample][["fq2"]]
    output:
        html="results/qc/fastqc/{sample}.2.html",
        zip="results/qc/fastqc/{sample}.2_fastqc.zip",
    log:
        "logs/fastqc/{sample}.2.log",
    threads: 8
    wrapper:
        "0.69.0/bio/fastqc"


# TODO include Kallisto
rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}.{read}_fastqc.zip", sample=get_samples(), read=[1,2]),
        expand(
            "results/species-diversity/{sample}/{sample}.uncleaned.kreport2", sample=get_samples()
        ),
        expand(
            "results/species-diversity-nonhuman/{sample}/{sample}.cleaned.kreport2", sample=get_samples()
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


# analysis of species diversity present BEFORE removing human contamination
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
        report="results/species-diversity/{sample}/{sample}.uncleaned.kreport2",
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


# plot Korna charts BEFORE removing human contamination
rule create_krona_chart:
    input:
        kraken_output="results/species-diversity/{sample}/{sample}.uncleaned.kreport2",
        taxonomy_database="resources/krona/",
    output:
        "results/species-diversity/{sample}/{sample}.html",
    log:
        "logs/krona/{sample}.log",
    conda:
        "../envs/kraken.yaml"
    shell:
        "ktImportTaxonomy -m 3 -t 5 -tax {input.taxonomy_database} -o {output} {input} 2> {log}"


# align trimmed reads against the human genome
rule align_against_human:
    input:
        "resources/genomes/human-genome.fna.gz",
        expand("results/trimmed/{{sample}}.{read}.fastq.gz", read=[1, 2]),
    output:
        "results/ordered-contigs-human/{sample}.bam",
    log:
        "logs/minimap2/{sample}.log",
    threads: 8
    resources:
        max_jobs = 1
    # TODO ask Johannes for usage if -x flag. Is "sr" ok?
    conda: 
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax sr -o {output} {input} 2> {log}"


# filter out human contamination
rule extract_unmapped:
    input:
        "results/ordered-contigs-human/{sample}.bam",
    output:
        fq1="results/ordered-contigs-nonhuman/{sample}.1.fastq.gz",
        fq2="results/ordered-contigs-nonhuman/{sample}.2.fastq.gz",
        sorted=temp("results/ordered-contigs-nonhuman/{sample}.bam"),
        unsorted=temp("results/ordered-contigs-nonhuman/{sample}.unsorted.bam"),
        t1=temp("results/ordered-contigs-nonhuman/{sample}.temp1.bam"),
        t2=temp("results/ordered-contigs-nonhuman/{sample}.temp2.bam"),
        t3=temp("results/ordered-contigs-nonhuman/{sample}.temp3.bam"),
    log:
        "logs/filter_human/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    threads: 8
    shell:
        """
        samtools view  -@ {threads} -u -f 4 -F 264 {input} > {output.t1}
        samtools view  -@ {threads} -u -f 8 -F 260 {input} > {output.t2}
        samtools view  -@ {threads} -u -f 12 -F 256 {input} > {output.t3}
        samtools sort  -@ {threads} -n -o {output.t1} {output.t1}
        samtools sort  -@ {threads} -n -o {output.t2} {output.t2}
        samtools sort  -@ {threads} -n -o {output.t3} {output.t3}
        samtools merge -@ {threads} -n {output.unsorted} {output.t1} {output.t2} {output.t3}
        samtools sort  -@ {threads} -n {output.unsorted} -o {output.sorted}
        samtools fastq -@ {threads} {output.sorted} -1 {output.fq1} -2 {output.fq2}
        """


# analysis of species diversity present AFTER removing human contamination
rule species_diversity_after:
    input:
        db="resources/minikraken-8GB",
        reads=expand(
            "results/ordered-contigs-nonhuman/{{sample}}.{read}.fastq.gz", read=[1, 2]
        ),
    output:
        kraken_output="results/species-diversity-nonhuman/{sample}/{sample}.kraken",
        report="results/species-diversity-nonhuman/{sample}/{sample}.cleaned.kreport2",
    log:
        "logs/kraken/{sample}_nonhuman.log",
    params:
    threads: 8
    conda:
        "../envs/kraken.yaml"
    shell:
        "(kraken2 --db {input.db} --threads {threads} --report {output.report} --gzip-compressed --paired {input.reads} > {output.kraken_output}) 2> {log}"


# plotting Krona charts AFTER removing human contamination
rule create_krona_chart_after:
    input:
        kraken_output="results/species-diversity-nonhuman/{sample}/{sample}.cleaned.kreport2",
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
