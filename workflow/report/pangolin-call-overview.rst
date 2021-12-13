Aggregated strain/lineage calls, as determined by `Pangolin <https://github.com/cov-lineages/pangolin>`_, throughout several stages of the workflow.

- Scaffolded Seq.: Contigs from De Novo assembly, ordered against the reference ({{ snakemake.config["virus-reference-genome"] }}) with `RaGOO <https://github.com/malonge/RaGOO>`_.
- Polished Seq.: Scaffolded sequence polished by applying variants with an allele frequency of 100% (called by `Varlociraptor <https://varlociraptor.github.io>`_ at FDR 5%).
- Masked Polished Seq.: Polished sequence masked (Positions in the reconstructed genome that are covered by less than {{ snakemake.config["quality-criteria"]["min-depth-with-PCR-duplicates"] }} reads must be masked with N. Informative positions (A, T, G, C) which are not supported by {{ snakemake.config["quality-criteria"]["min-allele"] }} of the aligned reads are masked by non-informative placeholders of the `IUPAC nucleotide code <https://www.bioinformatics.org/sms/iupac.html>`_.)
- Consensus Seq.: Only for Oxford Nanopore data. Consensus sequences of corrected reads and the reference genomes, as determined by `medaka <https://github.com/nanoporetech/medaka>`_.
- Masked Consensus Seq.: Only for Oxford Nanopore data. Masked consensus sequence. See above for details.
- Pseudo Seq.: Only for Illumina and Ion Torrent data. Sequence was generated based on the Wuhan Reference Genome ({{ snakemake.config["virus-reference-genome"] }}) and by applying variants with an allele frequency of 100% (called by `Varlociraptor <https://varlociraptor.github.io>`_ at FDR 5%).
