Strain calls for {{ snakemake.wildcards.sample }}. 
Horizontal axis shows the fraction in the sample, vertical axis shows the called strains.
Strains below {{ snakemake.params.min_fraction }} are summarized under the term "other".

Fraction estimation was performed with `Kallisto <https://pachterlab.github.io/kallisto>`_, taking the following predefined sequences as reference

{% for genome in snakemake.config["strain-calling"].get("genomes", []) %}
* {{ genome }}
{% endfor %}