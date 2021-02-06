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


rule multiqc:
    input:
        expand("results/qc/fastqc/{sample}_fastqc.zip", sample=get_samples()),
        expand("results/species-diversity/{sample}/{sample}.kreport2", sample=get_samples()),
    output:
        "results/qc/multiqc.html",
    params:
        "",  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc.log",
    wrapper:
        "0.69.0/bio/multiqc"


# Analysis of species diversity present
rule species_diversity:
    input:
        db = "resources/minikraken-8GB",
        reads = expand("results/trimmed/{{sample}}.{read}.fastq.gz", read=[1, 2])
    output:
        classified_reads = expand("results/species-diversity/{{sample}}/{{sample}}_{read}.classified.fasta", read=[1,2]),
        unclassified_reads = expand("results/species-diversity/{{sample}}/{{sample}}_{read}.unclassified.fasta", read=[1,2]), 
        kraken_output = "results/species-diversity/{sample}/{sample}.kraken",
        report = "results/species-diversity/{sample}/{sample}.kreport2"
    log:
        "logs/kraken/{sample}.log"
    params:
        classified = lambda w, output: "#".join(output.classified_reads[0].rsplit('_1', 1)),
        unclassified = lambda w, output: "#".join(output.unclassified_reads[0].rsplit('_1', 1))
    threads: 8
    conda:
        "../envs/kraken.yaml"
    shell:
        "(kraken2 --db {input.db} --threads {threads} --unclassified-out {params.unclassified} --classified-out {params.classified} --report {output.report} --gzip-compressed --paired {input.reads} > {output.kraken_output}) 2> {log}"

# Plot Korna charts
rule create_krona_chart:
    input:
        kraken_output = "results/species-diversity/{sample}/{sample}.kraken",
        taxonomy_database =  directory("resources/krona/")
    output:
        "results/species-diversity/{sample}/{sample}.html"
    log:
        "logs/krona/{sample}.log"
    conda:
        "../envs/kraken.yaml"
    shell:
        "ktImportTaxonomy -s 3 -t 4 -tax {input.taxonomy_database} -o {output} {input} 2> {log}"


# TODO Alexander and Thomas: add rules to detect contamination and perform QC
