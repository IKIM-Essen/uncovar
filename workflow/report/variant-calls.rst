Variant calls for {{ "all samples" if snakemake.wildcards.target == "all" else snakemake.wildcards.target }}, 
with false discovery rate (FDR) controlled at {{ snakemake.config["variant-calling"]["fdr"] }}, 
using the following filter expression: ``{{ snakemake.config["variant-calling"]["filters"][snakemake.wildcards.filter] }}``.

Candidate variants were obtained with `Freebayes <https://github.com/freebayes/freebayes>`_.
Variant calling and FDR control was performed with `Varlociraptor <https://varlociraptor.github.io>`_.