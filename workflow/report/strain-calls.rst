Strain calls for {{ snakemake.wildcards.sample }}. 
Horizontal axis shows the fraction in the sample, vertical axis shows the called strains.
Strains below {{ snakemake.params.min_fraction }} are summarized under the term "other".

Fraction estimation was performed with `Kallisto <https://pachterlab.github.io/kallisto>`_, taking all SARS-Cov-2 sequences from Genbank as reference.