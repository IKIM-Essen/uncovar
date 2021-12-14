# General settings
To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

# Sample sheet
The sample sheet contains all samples to be analyzed by UnCoVar.
## Auto filling

UnCoVar offers the possibility to automatically append samples to the sample sheet. To load your data into the workflow execute

    snakemake --cores all --use-conda update_sample

with the root of the UnCoVar as working directory. It is recommended to use the following structure to when adding data automatically:

    ├── archive
    ├── incoming
    └── snakemake-workflow-sars-cov2
        ├── data
        └── ...

However, this structure is not set in stone and can be adjusted via the `config/config.yaml` file under `data-handling`. Only the following path to the corresponding folders, relative to the directory of UnCoVar are needed:

- **incoming**: path of incoming data, which is moved to the data directory by the preprocessing script. Defaults to `../incoming/`.
- **data**: path to store data within the workflow. defaults to `data/`.
- **archive**: path to archive data from the results from the analysis to. Defaults to `../archive/`.

The incoming directory should contain paired end reads in (compressed) FASTQ format. UnCoVar automatically copies your data into the data directory and moves all files from incoming directory to the archive. After the analysis, all results are compressed and saved alongside the reads.

Moreover, the sample sheet is automatically updated with the new files. Please note, that only the part of the filename before the first '_' character is used as the sample name within the workflow.

## Manual filling

Of course, samples to be analyzed can also be added manually to the sample sheet. For each sample, the a new line in `config/pep/samples.csv` with the following content has to be defined:

- **sample_name**: name or identifier of sample
- **fq1**: path to read 1 in FASTQ format
- **fq2**: path to read 2 in FASTQ format
- **date**: sampling date of the sample
- **is_amplicon_data**: indicates whether the data was generated with a shotgun (0) or amplicon (1) sequencing
