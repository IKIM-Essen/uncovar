# get reference
rule hisat2_index:
    input:
        fasta="resources/genomes/main.fasta",
    output:
        directory("resources/genomes/indices/index_main"),
    params:
        extra="",
        prefix="resources/genomes/indices/index_main/",
    log:
        "logs/hisat2_index_main.log",
    threads: 8
    wrapper:
        "0.70.0/bio/hisat2/index"


# hisat2 rule
rule hisat2_align:
    input:
        reads=get_fastqs,
    output:
        "results/hisat2_mapped/{sample}.bam",
    log:
        "logs/hisat2_align_{sample}.log",
    params:
        extra="--dta-cufflinks",
        idx="resources/genomes/indices/index_main/",
    threads: 8
    wrapper:
        "0.70.0/bio/hisat2/align"


rule samtools_view:
    input:
        "results/hisat2_mapped/{sample}.bam",
    output:
        "results/hisat2_mapped/{sample}.sam",
    log:
        "logs/samtools/{sample}_view.log",
    params:
        "-h",  # optional params string
    wrapper:
        "0.70.0/bio/samtools/view"


rule samtools_flagstat:
    input:
        "results/hisat2_mapped/{sample}.sam",
    output:
        "results/hisat2_mapped/{sample}.sam.flagstat",
    log:
        "logs/samtools/{sample}_flagstat.log",
    wrapper:
        "0.70.0/bio/samtools/flagstat"


rule samtools_sort:
    input:
        "results/hisat2_mapped/{sample}.bam",
    output:
        "results/hisat2_mapped/{sample}.sorted.bam",
    log:
        "logs/samtools/{sample}_sort.log",
    params:
        extra="-m 16G",
        tmp_dir="/tmp/",
    # Samtools takes additional threads through its option -@
    threads: 8  # This value - 1 will be sent to -@.
    wrapper:
        "0.70.0/bio/samtools/sort"


rule generate_consensus:
    input:
        bam="results/hisat2_mapped/{sample}.sorted.bam",
        reference_genome="resources/genomes/main.fasta",
    output:
        "results/calls/{sample}_consensus.fq",
    log:
        "logs/bcftools_call/{sample}.log",
    conda:
        "../envs/sam+bcftools.yaml"
    shell:
        "(samtools mpileup -uf {input.reference_genome} {input.bam} |bcftools call -c |vcfutils.pl vcf2fq > {output}) &> {log}"


rule seqtk:
    input:
        "results/calls/{sample}_consensus.fq",
    output:
        "results/co_assembly_reads/{sample}.fasta",
    log:
        "logs/seqtk/{sample}.log",
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk seq -aQ64 -q20 -n N {input} > {output} 2> {log}"


rule call_variants:
    input:
        bam="results/hisat2_mapped/{sample}.sorted.bam",
        reference_genome="resources/genomes/main.fasta",
    output:
        "results/co_assembly_reads/{sample}_variants.vcf",
    log:
        "logs/bcftools_call/{sample}_vcf.log",
    conda:
        "../envs/sam+bcftools.yaml"
    shell:
        "(samtools mpileup -uf {input.reference_genome} {input.bam} |bcftools call -cv -Ob |bcftools view > {output}) &> {log}"


# evaluate the alignment statistics samtools flagstat
# samtools sort
# samtools faidx
# samtools mpileup
# seqtk
# samtools mpileup variant calling
# convert bcf to vcf
