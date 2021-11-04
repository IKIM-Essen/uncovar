# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    params:
        release=100,
    log:
        "logs/vep-plugins.log",
    wrapper:
        "0.69.0/bio/vep/plugins"


rule norm_bcfs:
    input:
        "results/{date}/calls/ref~main/{sample}.bcf",
    output:
        "results/{date}/calls/ref~main/normed_{sample}.bcf",
    log:
        "logs/{date}/norm-bcfs/{sample}.log",
    params:
        "-f /local/data/repos/snakemake-workflow-sars-cov2/resources/genomes/main.fasta -O b",
    wrapper:
        "0.79.0/bio/bcftools/norm"


rule annotate_variants:
    input:
        calls=get_bcf_for_annotation,
        plugins="resources/vep/plugins",
        fasta="resources/genomes/main.fasta",
        fai="resources/genomes/main.fasta.fai",
        gff="resources/annotation.gff.gz",
        csi="resources/annotation.gff.gz.tbi",
        problematic="resources/problematic-sites.vcf.gz",
        problematic_tbi="resources/problematic-sites.vcf.gz.tbi",
    output:
        calls="results/{date}/annotated-calls/ref~main/{normed_prefix}{sample}.bcf",
        stats="results/{date}/annotated-calls/ref~main/{normed_prefix}{sample}.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool"],
        extra=get_vep_args,
    log:
        "logs/{date}/vep/{normed_prefix}{sample}.log",
    wrapper:
        "0.72.0/bio/vep/annotate"
