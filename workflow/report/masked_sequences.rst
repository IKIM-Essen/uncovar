Masked, {{ snakemake.wildcards.reference }} sequence of sample {{ snakemake.wildcards.sample }}.

Positions in the reconstructed genome that are covered by less than {{ snakemake.wildcards.min_coverage }} reads are to be masked with N.
{% if snakemake.params.is_ont == False %}
Informative positions (A, T, G, C) supported by less than 90% of the aligned reads are masked by less or non-informative placeholders according to the `IUPAC nucleotide code <https://www.bioinformatics.org/sms/iupac.html>`_ usage (e.g. R, Y , N).
{% else %}
Since this sample was created using model-based basecall procedures (e.g. Oxford Nanopore), informative positions (A, T, G, C) which are not supported by 90% of the aligned reads are NOT masked.
{% endif %}