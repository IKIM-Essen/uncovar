# Understanding the Report

UnCoVar automatically generates a HTML-based report for a detailed insight into
the analysed patient or environmental samples.

## Section 1: Overview

### Report

An overview table of key information through several stages of the workflow.
The table consists of the following columns:

- Sample: The sample name as entered in the sample sheet
  (`config/pep/samples.csv`) before starting the workflow.
- Eukaryota to Unclassified: Classification of the sample contents on the
  domain level. Especially interesting for samples created via
  [shotgun sequencing](https://en.wikipedia.org/wiki/Shotgun_sequencing).
  with amplicon tiled preparation)
- Raw Reads: Number of reads in the unchanged sample.
- Trimmed Read: Number of reads after [adapter removal](https://www.ecseq.com/support/ngs/trimming-adapter-sequences-is-it-necessary).
- Filtered Reads: Number of reads after removing sequences that could belong
  to the human host from which the sample was taken.
- Best Quality: Indicates which of the assemblies performed produces the 
  highest quality sequence.
  filtering**and**primer clipping\*\*
- Largest Contig: The length of the largest sequences after
  [_de novo_ assembly](https://en.wikipedia.org/wiki/De_novo_sequence_assemblers)
- De Novo Sequence: The length of the _de novo_ assembled sequence after
  reference guided scaffolding. This sequence has less bias towards the
  SARS-CoV-2 reference genome.
  [RaGOO](https://github.com/malonge/RaGOO)
- Pseudo Sequence: The length of the pseudo assembled sequence. This sequence
  has a stronger bias towards the SARS-CoV-2 reference genome, as detected
  variants (on read level) are applied onto that reference. Only applicable
  for samples process on Illumina devices.
- Consensus Sequence: The length of the consensus sequences. This sequences
  also has a bias towards the SARS-CoV-2 reference genomes. It is created by
  using neural networks applied on a pileup of individual sequencing reads
  against the reference genome Only applicable for samples process on
  Oxford Nanopore devices.
  called variants with [Varlociraptor](https://varlociraptor.github.io)
  into the Wuhan-Hu-1 Reference
  [`NC_045512.2`](https://www.ncbi.nlm.nih.gov/nuccore/1798174254)
- Pango Lineage: Called lineage in [Pango nomenclature](https://cov-lineages.org/).
- WHO Label: Lineage name determined by the [WHO](https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/).
- VOC Mutations: Mutations that occur in Variants of Concern (VOC). These VOCs
  typically have: Increase in transmissibility or detrimental change in COVID-19
  epidemiology; OR Increase in virulence or change in clinical disease
  presentation; OR Decrease in effectiveness of public health and social
  measures or available diagnostics, vaccines, therapeutics.
  Other Mutations: Other Mutations that occur in the sample.
  of specific concern
  ([VOC](https://en.wikipedia.org/wiki/Variant_of_concern)) and other found variants

## Section 2: Variant Call Details

### VOC Similarity

Comparison of found mutations in the sample with mutations of Variants of
Concern. Show the top 5 Variants of Concern to the sample in terms of similarity.

- Mutations: All mutations found in the sample.
- Probability: A-posteriori probability that the mutations are present.
- Frequency: Variant allele frequency (VAF). The percentage of sequence reads
  observed matching a specific DNA variant divided by the overall coverage at
  that position and observed variant allele frequency (VAF) of each mutation
- Top 5 lineages: The most similar variants of concerns to the sample the
  [Jaccard index](https://en.wikipedia.org/wiki/Jaccard_index) is used as a
  measure of similarity.

## Section 3: Sequencing Details

## Section 4: Sequences

### Quality Overview

#### Filter Overview

#### Pangolin Call Overview

A table of lineage calls for all samples throughout several stages of the
workflow, as determined by [Pangolin](https://github.com/cov-lineages/pangolin).
These stages are:

- **Scaffolded** Sequences: Contigs from the
  [_de novo_ assembly](https://en.wikipedia.org/wiki/De_novo_sequence_assemblers),
  ordered against the virus reference genome (default: Wuhan-Hu-1 Reference
  [`NC_045512.2`](https://www.ncbi.nlm.nih.gov/nuccore/1798174254)) with
  [RaGOO](https://github.com/malonge/RaGOO).
- **Polished** Sequences: Scaffolded sequence polished by applying variants with
  an allele frequency of 100% (called by [Varlociraptor](https://varlociraptor.github.io)
  at FDR `5%` by default).
- **Masked Polished** Sequences: Polished Sequences masked due to two criteria.
  See [below](#quality-criteria) for a more in-depth explanation of these two criteria.
- **Consensus** Sequences: Using corrected reads and the virus reference
  genome (default: Wuhan-Hu-1 Reference
  [`NC_045512.2`](https://www.ncbi.nlm.nih.gov/nuccore/1798174254)) a
  [consensus sequence](https://en.wikipedia.org/wiki/Consensus_sequence)
  is generated using [medaka](https://github.com/nanoporetech/medaka).
  Consensus Sequences are only generated for Oxford Nanopore data.
- **Masked Consensus** Sequences: Consensus Sequence masked due to low position
  coverage. See [below](#quality-criteria) for a more in-depth explanation of
  the criteria. Masked consensus Sequences are only generated for Oxford
  Nanopore data.
- **Pseudo** Sequences: Sequences that are1 generated based on the virus reference
  genome (default: Wuhan-Hu-1 Reference
  [`NC_045512.2`](https://www.ncbi.nlm.nih.gov/nuccore/1798174254)) and by applying
  variants with an allele frequency of 100% (called by
  [Varlociraptor](https://varlociraptor.github.io) with a default FDR of `5%`).
  Masked Pseudo Sequences are not generated, as the quality criteria [below](#quality-criteria)
  are implicit considered when creating a "normal" Pseudo Sequence. Pseudo sequences
  are only generated for Illumina and Ion Torrent data.

##### Quality Criteria

1. Positions in the reconstructed genome that are covered by less than a certain
   amount (default: `20`) reads are be masked with N.
1. Informative positions (A, T, G, C) which are not supported by a certain
   percentage (default: `90%`) of the aligned reads are masked by non-informative
   placeholders of the [IUPAC nucleotide code](https://www.bioinformatics.org/sms/iupac.html).
   Not applicable for model-based basecall procedures (e.g. Oxford Nanopore).

## Section 5: Variant Call Files

## Section 6: High-Quality Genomes
