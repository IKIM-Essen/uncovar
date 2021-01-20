# TODO Johannes add rules to retrieve main assembly and annotation

rule get_genome:
    output:
        "refs/genome.fasta"
    log:
        "logs/get_genome.log"
    conda:
        "../envs/entrez.yaml"
    shell:
        "(esearch -db nucleotide -query 'NC_045512.2' |"
        "efetch -format fasta > {output}) 2> {log}"


rule get_genome_annotation:
    output:
        "refs/annotation.gff"
    log:
        "logs/get_annotation.log"
    conda:
        "../envs/curl.yaml"
    shell:
        "(curl -sSL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/"
        "GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz | "
        "zcat > {output}) 2> {log}"


rule get_problematic_sites:
    output:
        temp("refs/problematic-sites.vcf") # always retrieve the latest VCF
    conda:
        "../envs/curl.yaml"
    shell:
        "curl -sSL https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/"
        "master/problematic_sites_sarsCov2.vcf > {output} 2> {log}"

# TODO Alexander + Thomas add rules to retrieve strain sequences (I currently don't yet know from where)
