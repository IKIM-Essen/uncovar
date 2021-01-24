rule cat_genomes:
    input:
        get_strain_genomes,
    output:
        "resources/strain-genomes.fasta",
    log:
        "logs/cat-genomes.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "cat {input} > {output}"


# # NOT needed, if rule sourmash_compute_samples takes two fast qs
# rule cat_trimmed_samples:
#     input:
#         expand("results/trimmed/{{sample}}.{read}.fastq.gz", read=[1, 2])
#     output:
#         "results/trimmed/{sample}.both.fastq.gz"
#     shell:
#         "cat {input} > {output}"


rule sourmash_compute_genomes:
    input:
        "resources/strain-genomes.fasta",
    output:
        "resources/sourmash/genomes.sig",
    log:
        "logs/sourmash/sourmash-compute.log",
    threads: 4
    params:
        k="31",
        scaled="1000",
        extra="--singleton --track-abundance",  # compute signature for each sequence record individually
    wrapper:
        "0.70.0/bio/sourmash/compute"


rule sourmash_compute_samples:
    input:
        expand("results/trimmed/{{sample}}.{read}.fastq.gz", read=[1, 2]),
    output:
        "resources/sourmash/{sample}.sig",
    log:
        "logs/sourmash/sourmash-compute-{sample}.log",
    threads: 4
    params:
        k="31",
        scaled="1000",
        extra="--merge {sample} --track-abundance",
    wrapper:
        "0.70.0/bio/sourmash/compute"


# makes .sig from whole all single genomes .fasta
rule sourmash_compute:
    input:
        "resources/genomes/{accession}.fasta",
    output:
        "resources/genomes/{accession}.sig",
    log:
        "logs/sourmash/sourmash-compute-acc-{accession}.log",
    threads: 4
    params:
        k="31",
        scaled="1000",
        extra="",
    wrapper:
        "0.70.0/bio/sourmash/compute"


# db for smash search
rule smash_index:
    input:
        genomes=get_strain_signatures,
    output:
        "resources/sourmash/db.sbt.json",
    log:
        "logs/sourmash/index.log",
    conda:
        "../envs/sourmash.yaml"
    shell:
        "sourmash index -k 31 {output} {input.genomes}"


rule sourmash_search:
    input:
        read="resources/sourmash/{sample}.sig",
        db="resources/sourmash/db.sbt.json",
    output:
        "results/sourmash/search-{sample}.csv",
    log:
        "logs/sourmash/search-{sample}.log",
    conda:
        "../envs/sourmash.yaml"
    shell:
        "sourmash search {input.read} {input.db} -o {output} --threshold 0.001 --containment"


rule sourmash_gather:
    input:
        read="resources/sourmash/{sample}.sig",
        metagenome="resources/sourmash/genomes.sig",  # also vs single strain genome sig's?
    output:
        "results/sourmash/gather-{sample}.csv",
    log:
        "logs/sourmash/gather-{sample}.log",
    params:
        min_bp=config["strain-calling"]["min-bp"],
    conda:
        "../envs/sourmash.yaml"
    shell:
        "(sourmash gather -k 31 {input.read} {input.metagenome} -o {output} --threshold-bp {params.min_bp}) > {log} 2>&1"


# TODO
# 1.√at compute handover of two fastq files, dont use cat_trimmed_samples
# 2. the entrez rule get_genome also downloads partial reads of covid genomes (e.g. MW368461) or empty genome files (e.g. MW454604). Add rule to exiculde those rules "faulty" covid genomes
# 3.√Fix get_strain_accessions_from_txt - txt file must be in folder before staring the workflow -> questionable
# 3. clean up logging / add proper logging
# 4. if ok, then gather with folkers data
# 5. filter low abundance
# 6. parameter adjustment
# https://sourmash.readthedocs.io/en/latest/using-sourmash-a-guide.html#computing-signatures-for-read-files
