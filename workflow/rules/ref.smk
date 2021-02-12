checkpoint get_strain_accessions:
    output:
        "resources/strain-accessions.txt",
    log:
        "logs/get-accessions.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "curl -sSL https://www.ncbi.nlm.nih.gov/sars-cov-2/download-nuccore-ids > {output} 2> {log}"


rule get_genome:
    output:
        "resources/genomes/{accession}.fasta",
    params:
        accession=lambda w: "NC_045512.2" if w.accession == "main" else w.accession,
    log:
        "logs/genomes/get-genome/{accession}.log",
    conda:
        "../envs/entrez.yaml"
    resources:
        ncbi_api_requests=1,
    shell:
        "(esearch -db nucleotide -query '{params.accession}' | "
        "efetch -format fasta > {output}) 2> {log}"


rule get_genome_annotation:
    output:
        "resources/annotation.gff.gz",
    log:
        "logs/get-annotation.log",
    conda:
        "../envs/tabix.yaml"
    shell:
        # download, sort and bgzip gff (see https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html)
        "(curl -sSL ftp://ftp.ensemblgenomes.org/pub/viruses/gff3/sars_cov_2/Sars_cov_2.ASM985889v3.101.gff3.gz | "
        "zcat | grep -v '#' | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > {output}) 2> {log}"


rule get_problematic_sites:
    output:
        temp("resources/problematic-sites.vcf.gz"),  # always retrieve the latest VCF
    log:
        "logs/get-problematic-sites.log",
    conda:
        "../envs/tabix.yaml"
    shell:
        "curl -sSL https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/"
        "master/problematic_sites_sarsCov2.vcf | bgzip -c > {output} 2> {log}"
