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

# Tools, Frameworks and Packages used

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

