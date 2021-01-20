Variant calls for {{ "all samples" if snakemake.wildcards.target == "all" else snakemake.wildcards.target }}, 
with false discovery rate controlled at {{ snakemake.config["variant-calling"]["fdr"] }}.