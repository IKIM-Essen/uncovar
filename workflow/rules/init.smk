include: "ref.smk"

rule init:
    input: "resources/strain-accessions.txt",
           expand("resources/genomes/{accession}.fasta", accession=get_accessions()),
           "resources/annotation.gff.gz",
           "resources/protein_products.gff.gz",
           "resources/annotation_known_variants.gff.gz",
           "resources/minikraken-8GB/",
           "resources/krona/",
           "resources/genomes/human-genome.fna.gz",
           #temp("resources/gisaid/provision.json"),
           #expand("resources/genomes-renamed/{accession}.fasta", accession=get_accessions()),
           "data/",
           "../archive",
           "../incoming",


rule make_directories:
    output:
        data="data/",
        incoming="../incoming",
        archive="../archive"
    log:
        "log/make_directories.log"
    shell:
        "for dir in data/ ../archive/ ../incoming/; do if [ ! -d ""$dir"" ];"
        " then mkdir ""$dir""; fi done"