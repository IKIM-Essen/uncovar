# Snakemake workflow: SARS-Cov-2 strain and variant calling

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.27.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![GitHub actions status](https://github.com/koesterlab/snakemake-workflow-sars-cov-2/workflows/Tests/badge.svg?branch=master)](https://github.com/koesterlab/snakemake-workflow-sars-cov-2/actions?query=branch%3Amaster+workflow%3ATests)

This workflow calls the SARS-Cov-2 strain and variants for given NGS samples.

## Authors

* Alexander Thomas
* Thomas Battenfeld (@thomasbtf)
* Folker Meyer (@folker)
* Johannes Köster (@johanneskoester)

## Usage

### Step 1: Obtain a copy of this workflow

TODO upon publishing fill this with instructions on how to use the github template functionality.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files under `config`.

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N --resources ncbi_api_requests=1

using `$N` cores.
Non-local execution can be done via Snakemake's extensive cluster and cloud support, see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html).

### Step 4: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.zip