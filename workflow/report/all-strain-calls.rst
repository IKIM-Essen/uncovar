Strain calls for all samples. 
Horizontal axis shows the count of the respective strain, vertical axis shows the called strains.
{% if snakemake.wildcards.mode == "major" %}
Only the major strain in each sample is counted.
{% else %}
All strains above {{ snakemake.config["strain-calling"]["min-fraction"] * 100 }}% are counted.
{% endif %}

Fraction estimation was performed with `Kallisto <https://pachterlab.github.io/kallisto>`_, taking the following predefined sequences as reference

{% for genome in snakemake.config["strain-calling"].get("genomes", []) %}
* {{ genome }}
{% endfor %}