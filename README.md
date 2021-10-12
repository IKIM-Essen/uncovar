# UnCoVar: Snakemake workflow for SARS-Cov-2 strain and variant calling

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![GitHub actions status](https://github.com/koesterlab/snakemake-workflow-sars-cov2/workflows/Tests/badge.svg?branch=master)](https://github.com/koesterlab/snakemake-workflow-sars-cov2/actions?query=branch%3Amaster+workflow%3ATests)
[![Docker Repository on Quay](https://quay.io/repository/uncovar/uncovar/status "Docker Repository on Quay")](https://quay.io/repository/uncovar/uncovar)

![UnCoVar2](https://user-images.githubusercontent.com/77535027/133610563-d190e25c-504e-4953-92dd-f84a5b4a1191.png)

## Authors

* Alexander Thomas (@alethomas)
* Thomas Battenfeld (@thomasbtf)
* Felix Wiegand (@fxwiegand)
* Folker Meyer (@folker)
* Johannes Köster (@johanneskoester)

## Usage

### Step 1: Obtain a copy of this workflow

TODO upon publishing fill this with instructions on how to use the github template functionality.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files under `config`. Adjust `config.yaml` to configure the workflow execution, and `pep/samples.csv` to specify your sample setup.

#### Passing NGS samples
It is recommended to use the following structure to organize the data:

    ├── archive
    ├── incoming
    └── snakemake-workflow-sars-cov2
        ├── data
        └── ...

The incoming directory should contain paired end reads in a FASTQ format. It is recommended to work with compressed files (e.g. `sample-name.fastq.gz`).

To load your data into the workflow execute `python preprocessing/update_sample_sheet.py` with `snakemake-workflow-sars-cov2` as working directory.

The executed script automatically copies your data into the data directory and moves all files from incoming directory to the archive. 
Moreover, the sample sheet is automatically updated with the new files. Please note, that only the part of the filename before the first '_' character is used as the sample name for the workflow.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For Snakemake installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N --resources ncbi_api_requests=1

using `$N` cores.
Non-local execution can be done via Snakemake's extensive cluster and cloud support, see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html).

### Step 5: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.zip

## Tools, Frameworks and Packages used in UnCoVar

This project wouldn't be possible without several open source libraries:

| Tool           | Link                                              |
|----------------|---------------------------------------------------|
| ABySS          | www.doi.org/10.1101/gr.214346.116                 |
| Altair         | www.doi.org/10.21105/joss.01057                   |
| BAMClipper     | www.doi.org/10.1038/s41598-017-01703-6            |
| BCFtools       | www.doi.org/10.1093/gigascience/giab008           |
| BEDTools       | www.doi.org/10.1093/bioinformatics/btq033         |
| Biopython      | www.doi.org/10.1093/bioinformatics/btp163         |
| bwa            | www.doi.org/10.1093/bioinformatics/btp324         |
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
| mason          | www.http://publications.imp.fu-berlin.de/962      |
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
