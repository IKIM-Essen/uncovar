# The Report

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
- Raw Reads: Number of reads in the raw input of the sample.
- Trimmed Read: Number of reads after [adapter removal](https://www.ecseq.com/support/ngs/trimming-adapter-sequences-is-it-necessary).
- Filtered Reads: Number of reads after removing sequences that could belong
  to the human host from which the sample was taken.
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
- Consensus Sequence: The length of the consensus sequences. This sequence
  also has a bias towards the SARS-CoV-2 reference genomes. It is created by
  using neural networks applied on a pileup of individual sequencing reads
  against the reference genome Only applicable for samples process on
  Oxford Nanopore devices.
  called variants with [Varlociraptor](https://varlociraptor.github.io)
  into the Wuhan-Hu-1 Reference
  [`NC_045512.2`](https://www.ncbi.nlm.nih.gov/nuccore/1798174254)
- Best Quality: Indicates which of the assemblies performed produces the
  highest quality sequence and also indicates failed QC.
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

Estimation of VOC similarity based on each sample's SNV profile, compared to all
 VOCs from [covariants.org](covariants.org).

- Mutations: All possible SNVs found in the catalog with all SNVs of the most similar
  VOC showing first.
- Probability: A-posteriori probability that the observation of the mutation is true.
- Frequency: Variant allele frequency (VAF). The percentage of sequencing reads
  matching a specific SNV divided by the overall coverage at
  that position.
- The colmuns with the VOC names represent the Top 10 similar VOCs, based on the
  similarity of their SNV profile to the SNVs found in the sample.
  measure of similarity. To determine the degree of similarity between the sample
  and the VOCs, we performed the following scoring. Let $n$ be the total number
  of variants and $m$ be the number of VOCs. Let $X$ be a binary matrix that
  relates variants with VOCs, namely, $X_{i,j} = 1$ if and only if variant $i$
  is in  $VOC_j$, with $i = 1,\ldots,n$ and $j = 1,\ldots,m$. Let $\theta_i$
  be the latent allele frequency of variant $i$ in the given sample and $\hat\theta_i$
  be the maximum a posteriori estimate of $\theta_i$ as provided by Varlociraptor.
  Let $p_i = Pr(\theta_i > 0 \mid D)$ be the posterior probability that the variant
  $i$ is present in the sample (i.e., the probability that its latent allele
  frequency is greater than zero, given the data $D$). Then, the similarity of
  a given sample to $VOC_j$ can be calculated as the Jaccard-like similarity score:
  
$$\frac{\sum_{i=1}^n p_i \cdot \hat\theta_i \cdot X_{i,j} + (1-p_i) \cdot
 (1-X_{i,j})}{n}$$

### Rendered variant callings

Variant callings rendered with Oncoprint:

- with **high and moderate**/**low** impact
- on **ORF**/**Protein** level

Variant candidates are identified and called with Varlociraptor

## Section 3: Sequencing Details

### MultiQC report for overall quality control
- General Stats
- Sequence Counts
- Sequence Quality Histograms
- Per Sequence Quality Scores
- Per Base Sequence Content
- Per Sequence GC Content
- Per Base N Content
- Sequence Length Distribution
- Sequence Duplication Levels
- Overrepresented sequences by samples
- Top overrepresented sequences
- Adapter Content
- Status Checks
- Software Versions

### Coverage of Reference Genome

Plot that is showing how well the reference genome is covered by the reads of each
 analyzed sample visualizing the depth (DP) values from SAMtools depth.

## Section 4: Sequences

### Quality Overview

#### Filter Overview

A table comparing identity and number of N bases for all samples for the reconstructed
 genomes UnCoVar generates with its two main assembly methods (_de novo_
 assembly + scaffolding and consensus sequence based on called variants).

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

Downloadable `.vcf` files of the different stages for variant calling:

- with **high and moderate**/**low** impact
- on **ORF**/**Protein** level

## Section 6: High-Quality Genomes

- Multi-FASTA file, including the reconstructed genomes from samples, that passed
  the quality control
- `.csv` file, including additional submission data for each sample that passed
  the quality control
