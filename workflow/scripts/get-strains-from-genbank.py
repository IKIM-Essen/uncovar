import sys

sys.stderr = open(snakemake.log[0], "w")


def get_lineage_references():
    lineage_references = snakemake.params.lineage_references

    with open(snakemake.output[0], "w") as outfile:
        for _key, value in lineage_references.items():
            print("resources/genomes-renamed/{accession}.fasta".format(accession=value), file=outfile)

get_lineage_references()