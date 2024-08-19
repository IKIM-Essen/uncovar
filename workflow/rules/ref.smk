# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


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
        accession=(
            lambda w: config["virus-reference-genome"]
            if w.accession == "main"
            else w.accession
        ),
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


rule download_protein_products:
    output:
        temp("resources/protein_products.bed"),
    log:
        "logs/download_protein_products.log",
    conda:
        "../envs/ucsc.yaml"
    shell:
        "(bigBedToBed http://hgdownload.soe.ucsc.edu/gbdb/wuhCor1/uniprot/unipChainCov2.bb"
        " -chrom=NC_045512v2 -start=0 -end=29903 {output})"
        "2>{log}"


rule bed2gff:
    input:
        "resources/protein_products.bed",
    output:
        temp("resources/protein_products.gff"),
    log:
        "logs/bed2gff3.log",
    conda:
        "../envs/genometools.yaml"
    shell:
        "(cut -f1-12 {input} | sed -e 's/ /-/g' | sed -e 's/NC_045512v2/NC_045512.2/g'"
        " | gt bed_to_gff3 -featuretype gene -thicktype transcript -blocktype CDS -o {output} -force /dev/stdin )"
        "2> {log}"


rule filter_gff:
    input:
        "resources/protein_products.gff",
    output:
        temp("resources/protein_products.formatted.gff"),
    log:
        "logs/format_gff.log",
    conda:
        "../envs/tabix.yaml"
    shell:
        # download, sort and bgzip gff (see https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html)
        "cat {input} | grep -v '#' > {output} 2> {log}"


rule fix_gff:
    input:
        "resources/protein_products.formatted.gff",
    output:
        temp("resources/protein_products.fixed.gff"),
    log:
        "logs/fix_gff.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/fix-protein-gff.py"


rule format_fixed_gff:
    input:
        "resources/protein_products.fixed.gff",
    output:
        temp("resources/protein_products.fixed.formatted.gff"),
    log:
        "logs/format_gff.log",
    conda:
        "../envs/tabix.yaml"
    shell:
        # download, sort and bgzip gff (see https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html)
        "cat {input} | grep -v '#' | sort -k1,1 -k4,4n -k5,5n -k3,3n -t$'\t' > {output} 2> {log}"


rule compress_gff:
    input:
        "resources/protein_products.fixed.formatted.gff",
    output:
        "resources/protein_products.gff.gz",
    log:
        "logs/compress_gff.log",
    conda:
        "../envs/tabix.yaml"
    shell:
        "bgzip -c {input} > {output} 2> {log}"


rule get_genome_annotation_for_known_variants:
    output:
        "resources/annotation_known_variants.gff.gz",
    log:
        "logs/get-annotation_known_variants.log",
    conda:
        "../envs/tabix.yaml"
    shell:
        # download, sort and bgzip gff (see https://www.ensembl.org/info/docs/tools/vep/script/vep_custom.html)
        "(curl -sSL https://raw.githubusercontent.com/thomasbtf/nextclade/master/data/sars-cov-2/genemap.gff | "
        "cat | grep -v '#' | sort -k1,1 -k4,4n -k5,5n -t$'\t'  | bgzip -c > {output}) 2> {log}"


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
        "(mkdir {output} && curl -SL https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20220926.tar.gz | tar zxvf - -C {output}) 2> {log}"


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
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
        human_genome=config["human-genome-download-path"],
    log:
        "logs/get-human-genome.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "curl -SL -o {output} {params.human_genome} 2> {log}"


rule update_pangoLEARN:
    output:
        directory("results/{date}/pangolin/pangoLEARN"),
    log:
        "logs/{date}/pangolin/update.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {output} &&"
        " curl -L https://github.com/cov-lineages/pangoLEARN/archive/master.tar.gz |"
        " tar xvz --strip-components=1 -C {output})"
        " > {log} 2>&1"


rule update_lineages:
    output:
        directory("results/{date}/pangolin/lineages"),
    log:
        "logs/{date}/pangolin/update.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {output} &&"
        " curl -L https://github.com/cov-lineages/lineages/archive/master.tar.gz | "
        " tar xvz --strip-components=1 -C {output})"
        " > {log} 2>&1"


rule update_pangolin:
    output:
        "results/{date}/pangolin/update.log",
    log:
        "logs/{date}/pangolin/update.log",
    conda:
        "../envs/pangolin.yaml"
    shell:
        "pangolin --update-data > {log} 2>&1 && cp {log} {output}"
    


rule get_gisaid_provision:
    output:
        "resources/gisaid/provision.json",
    log:
        "logs/get_gisaid_provision.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(curl -L -u $GISAID_API_TOKEN https://www.epicov.org/epi3/3p/resseq02/export/provision.json.xz |"
        " xz -d -T0 > {output})"
        " > {log} 2>&1"


rule change_name_of_lineage_references:
    input:
        "resources/genomes/{accession}.fasta",
    output:
        "resources/genomes-renamed/{accession}.fasta",
    # get corresponding lineage of accession
    params:
        lineage=get_lineage_by_accession,
    log:
        "logs/change_name_of_lineage_references/{accession}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "sed -E 's/>(\S+)\\b/>{params.lineage}/;t' {input} > {output}"


checkpoint get_lineages_for_non_gisaid_based_calling:
    input:
        expand(
            "resources/genomes-renamed/{accession}.fasta",
            accession=config["strain-calling"]["lineage-references"].values(),
        ),
    output:
        "results/{date}/tables/predefinded-strain-genomes.txt",
    params:
        lineage_references=config["strain-calling"]["lineage-references"],
    log:
        "logs/{date}/get_lineages_for_non_gisaid_based_calling.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get-strains-from-genbank.py"
