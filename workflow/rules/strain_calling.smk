rule cat_covid_genomes:
    input:
        genomes = expand("resources/covid-genomes/{accession}.fasta", accession=get_strain_accessions_from_txt("resources/strain-accessions.txt"))
    output:
        "resources/covid-genomes.fasta"
    shell:
        "cat {input.genomes} > {output}"


rule cat_trimmed_samples:
    input:
        expand("results/trimmed/{{sample}}.{read}.fastq.gz", read=[1, 2])
    output:
        "results/trimmed/{sample}.both.fastq.gz"
    shell:
        "cat {input} > {output}"


rule sourmash_compute_covid_genomes:
    input:
        "resources/covid-genomes.fasta",
        # "resources/genomic.fna",
    output:
        "resources/sourmash/covid-genomes.sig",
    log:
        "logs/sourmash/sourmash-compute.log",
    threads: 2
    params:
        k="31",
        scaled="1000",
        extra="singleton",
    wrapper:
        "v0.69.0/bio/sourmash/compute"


rule sourmash_compute_samples:
    input:
        "results/trimmed/{sample}.both.fastq.gz"
    output:
        "resources/sourmash/{sample}.sig",
    log:
        "logs/sourmash/sourmash-compute-{sample}.log",
    threads: 2
    params:
        k="31",
        scaled="1000",
        extra="",
    wrapper:
        "v0.69.0/bio/sourmash/compute"


rule sourmash_gather:
    input:
        reads = "resources/sourmash/{sample}.sig",
        metagenome = "resources/sourmash/covid-genomes.sig"
    output:
        "results/sourmash/gather-{sample}.csv"
    log:
        "logs/sourmash/gather-{sample}.log"
    conda:
        "../envs/sourmash.yaml"
    shell:
        "(sourmash gather -k 31 {input.reads} {input.metagenome} -o {output}) 2> {log}"
