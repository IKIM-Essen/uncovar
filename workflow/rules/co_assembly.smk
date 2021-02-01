#get reference
rule hisat2_index:
    input:
        fasta = "resources/genomes/{sample}.fasta"
    output:
        directory("resources/genomes/indices/index_{sample}")
    params:
        extra = "",
        prefix = "resources/genomes/indices/index_{sample}/",
    log:
        "logs/hisat2_index_{sample}.log"
    threads: 8
    wrapper:
        "0.70.0/bio/hisat2/index"
        

# hisat2 rule
rule hisat2_align:
    input:
        reads=get_fastqs
    output:
        "results/hisat2_mapped/{sample}.bam"
    log:
        "logs/hisat2_align_{sample}.log"
    params:
      extra="--dta-cufflinks",
      idx="resources/genomes/indices/index_main",
    threads: 8
    wrapper:
      "0.70.0/bio/hisat2/align" 
# evaluate the alignment statistics samtools flagstat
# samtools sort
# samtools faidx
# samtools mpileup
# seqtk
# samtools mpileup variant calling
# convert bcf to vcf

