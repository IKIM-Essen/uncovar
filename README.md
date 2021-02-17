# Snakemake workflow: SARS-Cov-2 strain and variant calling

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.27.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![GitHub actions status](https://github.com/koesterlab/snakemake-workflow-sars-cov2/workflows/Tests/badge.svg?branch=master)](https://github.com/koesterlab/snakemake-workflow-sars-cov2/actions?query=branch%3Amaster+workflow%3ATests)

This workflow calls the SARS-Cov-2 strain and variants for given NGS samples.

## Authors

* Alexander Thomas (@alethomas)
* Thomas Battenfeld (@thomasbtf)
* Folker Meyer (@folker)
* Johannes Köster (@johanneskoester)

## Usage

### Step 1: Obtain a copy of this workflow

TODO upon publishing fill this with instructions on how to use the github template functionality.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files under `config`.

### Step 3: Passing data to the workflow

It is recommended to use the following structure to organize the data:

    ├── archive
    ├── incoming
    └── snakemake-workflow-sars-cov2
        ├── data
        └── ...

The incoming directory should contain paired end reads in a FASTQ format. It is recommended to work with compressed files (e.g. sample-name.fastq.gz).

To load your data into the workflow execute `python preprocessing/update_sample_sheet.py` from `snakemake-workflow-sars-cov2` as working directory.

The executed script automatically copies your data into the data directory of the workflow and moves all files from incoming directory to the archive. 
Moreover, the sample sheet is automatically updated with the new files. Please note, that only the part of the filename before the first '_' character is used as the sample name for the workflow.

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