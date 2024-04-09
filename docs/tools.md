# About

## The team

- Alexander Thomas
- Thomas Battenfeld
- Simon Magin
- Johannes Koester
- Folker Meyer

## Tools used in UnCoVar

| Stage | Step  | Tool Illumina | Tool Nanopore | SARS-CoV-2 specific |
|---|---|---|---|---|
| **Preprocessing** | primer clipping | [BAMClipper](https://github.com/tommyau/bamclipper) | [NoTrAmp](https://pypi.org/project/notramp/) | no |
|  | quality clipping | [fastp](https://github.com/OpenGene/fastp) | [fastp](https://github.com/OpenGene/fastp) | no |
|  | contamination removal | [Kraken2](https://github.com/DerrickWood/kraken2) | [Kraken2](https://github.com/DerrickWood/kraken2) | no |
|  | denoising | - | [Canu](https://github.com/marbl/canu), [Medaka](https://github.com/nanoporetech/medaka) | no |
| **Assembly** | assembly | [MEGAHIT](https://github.com/voutcn/megahit#), [metaSPAdes](https://github.com/ablab/spades) | [MEGAHIT](https://github.com/voutcn/megahit#), [metaSPAdes](https://github.com/ablab/spades) | no |
|  | scaffolding | [RaGOO](https://github.com/malonge/RaGOO) | [RaGOO](https://github.com/malonge/RaGOO) | no |
| **Analysis** | SNV calling | [freeBayes](https://github.com/freebayes/freebayes), [DELLY](https://github.com/dellytools/delly) | [Medaka](https://github.com/nanoporetech/medaka) variant, [Longshot](https://github.com/pjedge/longshot) | no |
|  | SNV validation | [Varlociraptor](https://github.com/varlociraptor/varlociraptor) | [Varlociraptor](https://github.com/varlociraptor/varlociraptor) | no |
| **Lineage call** | read based lineage assignment | [Kallisto](https://github.com/pachterlab/kallisto) | [Kallisto](https://github.com/pachterlab/kallisto) | no |
|  | lineage call | [Pangolin](https://github.com/cov-lineages/pangolin) | [Pangolin](https://github.com/cov-lineages/pangolin) | yes |

This project wouldn't be possible without several open source tools and libraries:

[Altair](www.doi.org/10.21105/joss.01057>)

[BCFtools](www.doi.org/10.1093/gigascience/giab008>)

[BEDTools](www.doi.org/10.1093/bioinformatics/btq033>)

[Biopython](www.doi.org/10.1093/bioinformatics/btp163>)

[bwa](www.doi.org/10.1093/bioinformatics/btp324>)

[Covariants](www.github.com/hodcroftlab/covariants>)

[delly](www.doi.org/10.1093/bioinformatics/bts378>)

[ensembl-vep](www.doi.org/10.1186/s13059-016-0974-4>)

[entrez-direct](www.ncbi.nlm.nih.gov/books/NBK179288>)

[FastQC](www.bioinformatics.babraham.ac.uk/projects/fastqc)

[fgbio](www.github.com/fulcrum-genomics/fgbio>)

[FreeBayes](www.arxiv.org/abs/1207.3907>)

[intervaltree](www.github.com/chaimleib/intervaltree>)

[Jupyter](www.jupyter.org>)

[kallisto](www.doi.org/10.1038/nbt.3519>)

[Krona](www.doi.org/10.1186/1471-2105-12-385>)

[mason](publications.imp.fu-berlin.de/962>)

[MEGAHIT](www.doi.org/10.1093/bioinformatics/btv033>)

[Minimap2](www.doi.org/10.1093/bioinformatics/bty191>)

[MultiQC](www.doi.org/10.1093/bioinformatics/btw354>)

[pandas](pandas.pydata.org>)

[Picard](broadinstitute.github.io/picard>)

[PySAM](www.doi.org/10.11578/dc.20190903.1>)

[QUAST](www.doi.org/10.1093/bioinformatics/btt086>)

[RaGOO](www.doi.org/10.1186/s13059-019-1829-6>)

[ruamel.yaml](www.sourceforge.net/projects/ruamel-yaml>)

[Rust-Bio-Tools](www.github.com/rust-bio/rust-bio-tools>)

[SAMtools](www.doi.org/10.1093/bioinformatics/btp352>)

[Snakemake](www.doi.org/10.12688/f1000research.29032.1>)

[sourmash](www.doi.org/10.21105/joss.00027>)

[SPAdes](www.doi.org/10.1089/cmb.2012.0021>)

[SVN](www.doi.org/10.1142/s0219720005001028>)

[Tabix](www.doi.org/10.1093/bioinformatics/btq671>)

[Trinity](www.doi.org/10.1038/nprot.2013.084>)

[Varlociraptor](www.doi.org/10.1186/s13059-020-01993-6>)

[Vega-Lite](www.doi.org/10.1109/TVCG.2016.2599030>)

[Velvet](www.doi.org/10.1101/gr.074492.107>)

[vembrane](www.github.com/vembrane/vembrane>)

