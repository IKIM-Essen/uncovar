Assembly of sample {{ snakemake.wildcards.sample }}. 
Reads were assembled with `Megahit <https://github.com/voutcn/megahit>`_, followed by reference based contig ordering and concatenation with `RaGOO <https://github.com/malonge/RaGOO>`_.
Then, assembly was polished by applying variants with an allele frequency of 100% (called by `Varlociraptor <https://varlociraptor.github.io>`_ at FDR 5%).