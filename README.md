# UnCoVar

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://github.com/IKIM-Essen/uncovar/assets/77535027/8e17c6fc-ff7a-4c25-afc9-7888036d693e" width="40%">
  <source media="(prefers-color-scheme: light)" srcset="https://github.com/IKIM-Essen/uncovar/assets/77535027/c99f5a94-749b-422e-b319-1e3700d40a8e" width="40%">
  <img alt="UnCoVar Logo dark/light">
</picture>

<h1>
Workflow for Transparent and Robust Virus Variant Calling, Genome Reconstruction
 and Lineage Assignment
</h1>

[![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![GitHub actions status](https://github.com/koesterlab/snakemake-workflow-sars-cov2/workflows/Tests/badge.svg?branch=master)](https://github.com/koesterlab/snakemake-workflow-sars-cov2/actions?query=branch%3Amaster+workflow%3ATests)
[![Docker Repository on Quay](https://quay.io/repository/uncovar/uncovar/status)](https://quay.io/repository/uncovar/uncovar)

## Usage

The usage is described in detail in UnCoVar's [GitHub pages](https://ikim-essen.github.io/uncovar/)

### Step 1: Install Snakemake and Snakedeploy

Snakemake and Snakedeploy are best installed via the [Mamba package manager](https://github.com/mamba-org/mamba)
 (a drop-in replacement for conda). If you have neither Conda nor Mamba, it can
  be installed via [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge).
  For other options see [here](https://github.com/mamba-org/mamba).

Given that Mamba is installed, run

```sh
    mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
```

to install both Snakemake and Snakedeploy in an isolated environment. For all
 following commands ensure that this environment is activated via

```sh
    conda activate snakemake
```

### Step 2: Clone or Deploy workflow

First, create an appropriate project working directory on your system and enter it:

```sh
    WORKDIR=path/to/project-workdir
    mkdir -p ${WORKDIR}
    cd ${WORKDIR}
```

In all following steps, we will assume that you are inside of that directory.
Second, run

Given that Snakemake is installed and you want to clone the full workflow you can
 do it as follows:

```sh
    git clone https://github.com/IKIM-Essen/uncovar
```

Given that Snakemake and Snakedeploy are installed and available (see Step 1),
 the workflow can be deployed as follows:

```sh
    snakedeploy deploy-workflow https://github.com/IKIM-Essen/uncovar . --tag v0.16.0
```

Snakedeploy will create two folders `workflow` and `config`. The former contains
 the deployment of the UnCoVar workflow as a
  [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#using-and-combining-pre-exising-workflows),
  the latter contains configuration files which will be modified in the next step
  in order to configure the workflow to your needs. Later, when executing the workflow,
  Snakemake will automatically find the main Snakefile in the workflow subfolder.

### Step 3: Configure workflow

### General settings

#### Config file

To configure this workflow, modify `config/config.yaml` according to your
 needs, following the explanations provided in the file. It is especially recommended
 to provide a `BED` file with primer coordinates, when the analyzed reads derive
 from amplicon-tiled sequencing, so the primers are trimmed appropriately.

#### Sample sheet

The sample sheet contains all samples to be analyzed by UnCoVar.

#### Auto filling

UnCoVar offers the possibility to automatically append paired-end sequenced
 samples to the sample sheet. To load your data into the workflow execute

```sh
    snakemake --cores all --use-conda update_sample
```

with the root of the UnCoVar as working directory. It is recommended to use
the following structure to when adding data automatically:

```text
    ├── archive
    ├── incoming
    └── snakemake-workflow-sars-cov2
        ├── data
            └── 2023-12-24
```

However, this structure is not set in stone and can be adjusted via the
`config/config.yaml` file under `data-handling`. Only the following path to the
corresponding folders, relative to the directory of UnCoVar are needed:

- **incoming**: path of incoming data, which is moved to the data directory by
  the preprocessing script. Defaults to `../incoming/`.
- **data**: path to store data within the workflow. defaults to `data/`. It is
 recommend using subfolders named properly (e.g. with date)
- **archive**: path to archive data from the results from the analysis to.
  Defaults to `../archive/`.

The incoming directory should contain paired end reads in (compressed) FASTQ
format. UnCoVar automatically copies your data into the data directory and moves
all files from incoming directory to the archive. After the analysis, all results
are compressed and saved alongside the reads.

Moreover, the sample sheet is automatically updated with the new files. Please
 note, that only the part of the filename before the first '\_' character is used
 as the sample name within the workflow. Technology and amplicon flag (**is_amplicon_data**)
 have to be revisited manually

#### Manual filling

Of course, samples to be analyzed can also be added manually to the sample sheet.
For each sample, the a new line in `config/pep/samples.csv` with the following
content has to be defined:

- **sample_name**: name or identifier of sample
- **fq1**: path to read 1 in FASTQ format
- **fq2**: path to read 2 in FASTQ format (if paired end sequencing)
- **date**: sampling date of the sample
- **is_amplicon_data**: indicates whether the data was generated with a
  shotgun (0) or amplicon (1) sequencing
- **technology**: indicates the sequencing technology used to generate
  the samples (illumina, ont, ion)
- **include_in_high_genome_summary**: indicates if sample should be included in the submission files (1) or not (0)

### Step 4: Run workflow

Given that the workflow has been properly deployed and configured, it can be executed as follows.

Fow running the workflow while deploying any necessary software via conda (using
 the Mamba package manager by default), run Snakemake with

```sh
    snakemake --cores all --use-conda
```

Snakemake will automatically detect the main Snakefile in the workflow subfolder
 and execute the workflow module that has been defined by the deployment in step 2.

For further options, e.g. for cluster and cloud execution, see the docs.


This workflow is written with Snakemake and details and tools are described in the
[Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=IKIM-Essen/uncovar).

If you use this workflow in a paper, don't forget to give credits to the
authors by citing the URL of this repository and its DOI (see above).
