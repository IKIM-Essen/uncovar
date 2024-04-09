# Advanced Configuration

The config file, found under `config/config.yaml` can be used to adapt your analysis.

## Execution Mode

```yaml
# execution mode. Can be either "patient" or "environment"
mode: environment
```

Defines the execution mode of UnCoVar.

When the mode is set to `patient`, the sample is assumed come be from a single
host organism and contains only one strain of SARS-CoV-2. The parts of the
workflow for reconstructing the SARS-CoV-2 strain genome are activated.

If the mode is set to `environment`, the sample is assumed to be from the
environment (e.g. wastewater) and to contain different SARS-CoV-2 strains.
The parts of the workflow responsible for creating and analysing individual
genomes (e.g. assembly, lineage calling via Pangolin) are disabled.

## Sending lab number

UnCoVar automatically generates a multi-Fasta file and a corresponding `.csv` for
 all samples with a `1` flag for `inlcude_in_high_genome_summary` in the sample sheet,
 that match the given `quality-criteria` (see below). The reporting format and the
 quality criteria are inspired by the [requirements for SARS-CoV-2 genome submission](https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/DESH/Qualitaetskriterien.pdf?__blob=publicationFile)
 to the [Robert-Koch-Institute, Germany](https://www.rki.de/DE/Content/InfAZ/N/Neuartiges_Coronavirus/nCoV.html).
 The sending lab number will be included in the `.csv` file

## Data handling

With the root of the UnCoVar workflow as working directory, we recommended to
 use the following folder structure:

```text
├── archive
├── incoming
└── uncovar
    └── data
        └── 2023-12-24
```

The structure can be adjusted to via the config under `data-handling`:

```yaml
data-handling:
  # flag for using the following data-handling structure
  # True: data-handling structure is used as shown below
  # False: only the sample sheet needs to be updated (manually)
  use-data-handling: True
  # flag for archiving data
  # True: data is archived in path defined below
  # False: data is not archived
  archive-data: False
  # path of incoming data, which is moved to the
  # data directory by the preprocessing script
  incoming: ../incoming/
  # path to store data within the workflow
  data: data/
  # path to archive data from incoming and
  # the results from the latest run to
  archive: ../archive/
```

## Quality criteria

The quality criteria can be adjusted to your individual needs. By default they match
 the quality criteria needed for submitting to the RKI (see **Sending lab number**
 above)

```yaml
quality-criteria:
  illumina:
    # minimal length of acceptable reads
    min-length-reads: 30
    # average quality of acceptable reads (PHRED)
    min-PHRED: 20
  ont:
    # minimal length of acceptable reads
    min-length-reads: 200
    # average quality of acceptable reads (PHRED)
    min-PHRED: 10
  # identity to virus reference genome (see-above) of reconstructed sequence
  min-identity: 0.9
  # share N in the reconstructed sequence
  max-n: 0.05
  # minimum local sequencing depth without filtering of PCR duplicates
  min-depth-with-PCR-duplicates: 20
  # minimum local sequencing depth after filtering PCR duplicates
  min-depth-without-PCR-duplicates: 10
  # minimum informative allele frequency
  min-allele: 0.9
```

## Preprocessing

Here different preprocessing can be adjustet. Per default the standard Illumina adapters
 are trimmed. For samples prepared with an amplicon sequencing approach, you can
 define the path to the primer file in `.bed` format. If you are processing Nanopore
 samples, you can also define the primer version via changing the number.

The default primer file is a bed file from the [ARTIC network](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3>).
However, the primers for clipping can be customized. First, the custom primers must
be saved in bed format. Next, the path to this file must be changed in the config.
Go to the config folder and open config.yaml. In the "preprocessing" subcategory,
change the path after "amplicon-primers" to the path where your primer file
can be found.

```yaml
preprocessing:
  # only for *non* Oxford Nanopore data. Adapters to trim.
  # see: https://www.nimagen.com/shop/products/rc-cov096/easyseq-sars-cov-2-novel-coronavirus-whole-genome-sequencing-kit
  kit-adapters: "--adapter_sequence GCGAATTTCGACGATCGTTGCATTAACTCGCGAA --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
  # only for Oxford Nanopore data.
  # ARTIC primer version to clip from reads. See
  # https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V4
  # for more information
  artic-primer-version: 4
  # path to amplicon primers in bed format for hard-clipping on paired end files (illumina) or url to file that should be downloaded
  amplicon-primers: "resources/SARS-CoV-2-artic-v4_1.primer.bed"
  # GenBank accession of reference sequence of the amplicon primers
  amplicon-reference: "MN908947"
```

## Assembly

In this section you define which assembler you want to use for the genome reconstruction.
 UnCoVar uses MEGAHIT and metaSPAdes by default, as those achieved the best results
 in a benchmarking comparison. The assembly options can be changed independently.

There are several other options available:

- megahit-std
- megahit-meta-large
- megahit-meta-sensitive
- trinity
- velvet
- metaspades
- coronaspades
- spades
- rnaviralspades
