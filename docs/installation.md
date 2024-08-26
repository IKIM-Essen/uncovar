# Installation

## Install Snakemake and Snakedeploy

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

## Clone or deploy workflow

First, create an appropriate project working directory on your system and enter it:

```sh
    WORKDIR=path/to/project-workdir
    mkdir -p ${WORKDIR}
    cd ${WORKDIR}
```

In all following steps, we will assume that you are inside of that directory.

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

## Sample and workflow preparation

UnCoVar is now ready to be used and will download all neccesary resources automatically.
 To start your first analysis, you will need input files in FASTQ-format, that are
 defined in the sample sheet in: `config/pep/samples.csv`. If you are using
 read files generated with an amplicon protocol, you will also need to provide
 the amplicon primer coordinates in `BED` format, so the primers are trimmed appropriately.

With the root of the UnCoVar workflow as working directory, we recommended to
 use the following folder structure:

```text
    ├── archive
    ├── incoming
    └── uncovar
        └── data
            └── 2023-12-24
```

The structure can be adjusted via the config file: `config/config.yaml` under
 `data-handling`:

- **incoming**: path of incoming data, which is moved to the data directory by
  the preprocessing script. Defaults to `../incoming/`.
- **data**: path to store data within the workflow. defaults to `data/`. It is
 recommend using subfolders named properly (e.g. with date)
- **archive**: path to archive data from the results from the analysis to.
  Defaults to `../archive/`.

### ONT sample preprocessing
Automates the extraction of sequencing reads from sample-specific folders, merges the reads from each 
  barcode into a single FASTQ file, and then renames the files for easy identification and downstream analysis.

- **barcode-rename.csv**: A CSV file containing the barcode sequences and their corresponding sample 
  identifiers.
- **barcode_dir**: The directory path where the barcode-specific folders are located.
- **output_dir**: The directory path where the renamed FASTQ files will be saved.

To run the tool, use the following command:
```sh
  snakemake --config barcode_dir=path/to/barcode/folders output_dir=data/date/ --cores all --use-conda barcode_rename
```

### Sample sheet

The sample sheet contains all samples to be analyzed by UnCoVar. UnCoVar offers
 the possibility to automatically append paired-end sequenced
 samples to the sample sheet. Single-end sequenced samples have to be added manually
 (see **Manual filling** below).

### Auto filling

To UnCoVars automated sample sheet filling, place your raw and compressed
 FASTQ-files into the `../incoming` folder and run:

```sh
    snakemake --cores all --use-conda update_sample
```

The files are appended to the sample sheet (auto check for duplications) and are
 moved into a new folder with the current date (`YYYY-MM-DD`). `sample_name`
 is extracted from each filename and all characters before the first \_ are used.
 Make sure to include the correct flags for amplicon-generated files
 (`is_amplicon_data`) and sequencing technology (`technology`). For details on
 the flags for each field see below.

### Manual filling

Of course, samples to be analyzed can also be added manually to the sample sheet.
 For each sample a new line in `config/pep/samples.csv` with the following
 information has to be defined:

| sample_name | fq1 | fq2 | date | is_amplicon_data | technology | include_in_high_genome_summary |
| --- | --- | --- | --- | --- | --- | --- |
| example-1 | PATH/TO/fq1 | PATH/TO/fq2 | 2024-01-01 | 1 | illumina | 1 |
| example-2 | PATH/TO/fq | | 2024-01-01 | 1 | ont | 1 |

- **sample_name**: name or identifier of sample
- **fq1**: path to read 1 in compressed FASTQ format
- **fq2**: path to read 2 in conpressed FASTQ format (if paired-end sequencing)
- **date**: sampling date of the sample
- **is_amplicon_data**: indicates whether the data was generated with a
  shotgun (0) or amplicon (1) sequencing
- **technology**: indicates the sequencing technology used to generate
  the samples (illumina, ont, ion)
- **include_in_high_genome_summary**: indicates if sample should be included in
 the submission files (1) or not (0)

## Run the workflow

Given that the workflow has been properly deployed and configured, run Snakemake
 with:

```sh
    snakemake --cores all --use-conda
```

Snakemake will automatically detect the main Snakefile in the workflow subfolder
 and execute the workflow module that has been defined by the deployment.

This workflow is written with Snakemake and details and tools are described in the
[Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=IKIM-Essen/uncovar).

If you use this workflow in your work, don't forget to give credits to the
authors by citing the URL of this repository and its DOI (see above).

----------------
## General settings

### Config file

To configure this workflow, modify `config/config.yaml` according to your
 needs, following the explanations provided in the file. It is especially recommended
 to provide a `BED` file with primer coordinates, when the analyzed reads derive
 from amplicon-tiled sequencing, so the primers are trimmed appropriately.

The incoming directory should contain paired end reads in (compressed) FASTQ
format. UnCoVar automatically copies your data into the data directory and moves
all files from incoming directory to the archive. After the analysis, all results
are compressed and saved alongside the reads.
