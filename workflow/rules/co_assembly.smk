#get reference
rule hisat2_index:
    input:
        fasta = "resources/genomes/main.fasta"
    output:
        directory("resources/genomes/indices/index_main")
    params:
        prefix = "index_main/"
    log:
        "logs/hisat2_index_main.log"
    threads: 8
    wrapper:
        "0.70.0/bio/hisat2/index"
        


# hisat2 rule
# rule hisat:
#     input:
        
#         reads=get_fastqs
#     output:
#         "results/mapped/{sample}.bam"
#     log:
#         "logs/hisat2_align_{sample}.log"
#     params:
#       extra="-p 16 -dta-cufflinks",
#       idx="index/",
#     threads: 8
#     wrapper:
#       "0.70.0/bio/hisat2/align" 
# evaluate tghe alignment statistics samtools flagstat
# samtools sort
# samtools faidx
# samtools mpileup
# seqtk
# samtools mpileup variant calling
# convert bcf to vcf

