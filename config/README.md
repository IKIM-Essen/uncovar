# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet

* Manual filling: Add samples to `config/pep/samples.csv`. For each sample, the columns `sample_name`, `fq1`, `fq2`,`run_id`, `is_amplicon_data` have to be defined.
* Automatic filling: It is recommended to use the following structure to when adding data via the provided preprocessing script:

    ├── archive
    ├── incoming
    └── snakemake-workflow-sars-cov2
        ├── data
        └── ...

The incoming directory should contain paired end reads in (compressed) FASTQ format. 
To load your data into the workflow execute `python preprocessing/update_sample_sheet.py` with the root of the workflow as working directory.

The executed script automatically copies your data into the data directory and moves all files from incoming directory to the archive. Moreover, the sample sheet is automatically updated with the new files. Please note, that only the part of the filename before the first '_' character is used as the sample name for the workflow.
