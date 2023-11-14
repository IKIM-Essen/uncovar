# UnCoVar


<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/IKIM-Essen/uncovar/assets/77535027/c99f5a94-749b-422e-b319-1e3700d40a8e">
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/IKIM-Essen/uncovar/assets/77535027/8e17c6fc-ff7a-4c25-afc9-7888036d693e">
  <img alt="UnCoVar Logo dark/light">
</picture>


## SARS-CoV-2 Variant Calling and Lineage Assignment

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![GitHub actions status](https://github.com/koesterlab/snakemake-workflow-sars-cov2/workflows/Tests/badge.svg?branch=master)](https://github.com/koesterlab/snakemake-workflow-sars-cov2/actions?query=branch%3Amaster+workflow%3ATests)
[![Docker Repository on Quay](https://quay.io/repository/uncovar/uncovar/status)](https://quay.io/repository/uncovar/uncovar)

A reproducible and scalable workflow for transparent and robust SARS-CoV-2
variant calling and lineage assignment with comprehensive reporting.

## Usage

This workflow is written with snakemake and its usage is described in the
[Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=IKIM-Essen/uncovar).

If you use this workflow in a paper, don't forget to give credits to the
authors by citing the URL of this repository and its DOI (see above).

## Tools, Frameworks and Packages used

This project wouldn't be possible without several open source libraries:

| Tool           | Link                                              |
| -------------- | ------------------------------------------------- |
| ABySS          | www.doi.org/10.1101/gr.214346.116                 |
| Altair         | www.doi.org/10.21105/joss.01057                   |
| BAMClipper     | www.doi.org/10.1038/s41598-017-01703-6            |
| BCFtools       | www.doi.org/10.1093/gigascience/giab008           |
| BEDTools       | www.doi.org/10.1093/bioinformatics/btq033         |
| Biopython      | www.doi.org/10.1093/bioinformatics/btp163         |
| bwa            | www.doi.org/10.1093/bioinformatics/btp324         |
| Covariants     | www.github.com/hodcroftlab/covariants             |
| delly          | www.doi.org/10.1093/bioinformatics/bts378         |
| ensembl-vep    | www.doi.org/10.1186/s13059-016-0974-4             |
| entrez-direct  | www.ncbi.nlm.nih.gov/books/NBK179288              |
| fastp          | www.doi.org/10.1093/bioinformatics/bty560         |
| FastQC         | www.bioinformatics.babraham.ac.uk/projects/fastqc |
| fgbio          | www.github.com/fulcrum-genomics/fgbio             |
| FreeBayes      | www.arxiv.org/abs/1207.3907                       |
| intervaltree   | www.github.com/chaimleib/intervaltree             |
| Jupyter        | www.jupyter.org                                   |
| kallisto       | www.doi.org/10.1038/nbt.3519                      |
| Kraken2        | www.doi.org/10.1186/s13059-019-1891-0             |
| Krona          | www.doi.org/10.1186/1471-2105-12-385              |
| mason          | www.<http://publications.imp.fu-berlin.de/962>    |
| MEGAHIT        | www.doi.org/10.1093/bioinformatics/btv033         |
| Minimap2       | www.doi.org/10.1093/bioinformatics/bty191         |
| MultiQC        | www.doi.org/10.1093/bioinformatics/btw354         |
| pandas         | pandas.pydata.org                                 |
| Picard         | broadinstitute.github.io/picard                   |
| PySAM          | www.doi.org/10.11578/dc.20190903.1                |
| QUAST          | www.doi.org/10.1093/bioinformatics/btt086         |
| RaGOO          | www.doi.org/10.1186/s13059-019-1829-6             |
| ruamel.yaml    | www.sourceforge.net/projects/ruamel-yaml          |
| Rust-Bio-Tools | www.github.com/rust-bio/rust-bio-tools            |
| SAMtools       | www.doi.org/10.1093/bioinformatics/btp352         |
| Snakemake      | www.doi.org/10.12688/f1000research.29032.1        |
| sourmash       | www.doi.org/10.21105/joss.00027                   |
| SPAdes         | www.doi.org/10.1089/cmb.2012.0021                 |
| SVN            | www.doi.org/10.1142/s0219720005001028             |
| Tabix          | www.doi.org/10.1093/bioinformatics/btq671         |
| Trinity        | www.doi.org/10.1038/nprot.2013.084                |
| Varlociraptor  | www.doi.org/10.1186/s13059-020-01993-6            |
| Vega-Lite      | www.doi.org/10.1109/TVCG.2016.2599030             |
| Velvet         | www.doi.org/10.1101/gr.074492.107                 |
| vembrane       | www.github.com/vembrane/vembrane                  |
