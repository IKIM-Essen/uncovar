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
           expand("resources/genomes-renamed/{accession}.fasta", accession=config["strain-calling"]["lineage-references"].values()),
           ancient("data/"),
           ancient("../archive"),
           ancient("../incoming"),

rule make_directory_data:
    output:
        data=directory("data"),
    log:
        "logs/make_directory_data.log"
    shell:
        "for dir in {output}; do if [ ! -d ""$dir"" ];"
        " then mkdir ""$dir""; fi done"

rule make_directory_incoming:
    output:
        incoming=directory("../incoming/"),
    log:
        "logs/make_directory_incoming.log"
    shell:
        "for dir in {output}; do if [ ! -d ""$dir"" ];"
        " then mkdir ""$dir""; fi done"

rule make_directory_archive:
    output:
        archive=directory("../archive"),
    log:
        "logs/make_directory_archive.log"
    shell:
        "for dir in {output}; do if [ ! -d ""$dir"" ];"
        " then mkdir ""$dir""; fi done"
    