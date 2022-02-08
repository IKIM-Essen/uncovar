rule minimap2_bam_sanger:
    input:
        target="resources/genomes/main.fasta",
        query=get_sanger_files_for_sample,
    output:
        "results/benchmarking/sanger/aligned/{sample}.bam",
    log:
        "logs/minimap2/{sample}.log",
    params:
        extra="",
        sorting="coordinate",
        sort_extra="",
    threads: 3
    wrapper:
        "v1.0.0/bio/minimap2/aligner"


rule bamtobed:
    input:
        "results/benchmarking/sanger/aligned/{sample}.bam",
    output:
        "results/benchmarking/sanger/aligned/{sample}.bed",
    log:
        "logs/bamtobed/{sample}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "bamToBed -i {input} > {output} 2> {log}"


rule freebayes_sanger:
    input:
        # optional BED file specifying chromosomal regions on which freebayes
        # should run, e.g. all regions that show coverage
        #regions="path/to/region-file.bed"
        ref="resources/genomes/main.fasta",
        fai="resources/genomes/main.fasta.fai",
        samples="results/benchmarking/sanger/aligned/{sample}.bam",
        indexes="results/benchmarking/sanger/aligned/{sample}.bam.bai",
    output:
        "results/benchmarking/sanger/variant-calls/{sample}.vcf",  # either .vcf or .bcf
    log:
        "logs/freebayes/{sample}.log",
    params:
        # TODO: GT 1 or 1/1
        extra="--min-alternate-count 1 --ploidy 2 --pooled-continuous",  # optional parameters
        chunksize=100000,  # reference genome chunk size for parallelization (default: 100000)
        normalize=False,  # optional flag to use bcftools norm to normalize indels (Valid params are -a, -f, -m, -D or -d)
    threads: 2
    wrapper:
        "v1.0.0/bio/freebayes"
