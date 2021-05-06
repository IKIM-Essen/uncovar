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
        "((esearch -db nucleotide -query '{params.accession}' | "
        "efetch -format fasta > {output}) && [ -s {output} ]) 2> {log}"


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


rule get_genome_db_for_kraken:
    output:
        directory("resources/minikraken-8GB"),
    log:
        "logs/get-kraken-db.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "mkdir {output} && curl -SL ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/minikraken_8GB_202003.tgz | tar zxvf - -C {output} --strip 1 2> {log}"


rule get_taxonomie_db_for_krona:
    output:
        directory("resources/krona/"),
    log:
        "logs/get-krona-db.log",
    conda:
        "../envs/kraken.yaml"
    shell:
        "ktUpdateTaxonomy.sh {output} 2> {log}"


rule get_human_genome:
    output:
        "resources/genomes/human-genome.fna.gz",
    log:
        "logs/get-human-genome.log",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    conda:
        "../envs/unix.yaml"
    shell:
        "curl -SL -o {output} ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz 2> {log}"


rule update_pangoLEARN:
    output:
        directory("results/{date}/pangolin/pangoLEARN"),
    log:
        "logs/{date}/pangolin/update.log",
    shell:
        "(mkdir -p {output} && wget -qO- https://github.com/cov-lineages/pangoLEARN/archive/master.tar.gz | tar xvz --strip-components=1 -C {output})> {log} 2>&1"


rule update_lineages:
    output:
        directory("results/{date}/pangolin/lineages"),
    log:
        "logs/{date}/pangolin/update.log",
    shell:
        "(mkdir -p {output} && wget -qO- https://github.com/cov-lineages/lineages/archive/master.tar.gz | tar xvz --strip-components=1 -C {output})> {log} 2>&1"


rule get_gisaid_provision:
    output:
        temp("resources/gisaid/provision.json"),
    log:
        "logs/get_gisaid_provision.log",
    shell:
        "(curl -u $GISAID_API_TOKEN  https://www.epicov.org/epi3/3p/resseq02/export/provision.json.xz | xz -d -T0 > {output})> {log} 2>&1"
