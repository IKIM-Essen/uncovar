# Understanding the Report

## Section 1: Overview

## Section 2: Variant Call Details

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
