rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    params:
        release=100,
    log:
        "logs/vep-plugins.log",
    wrapper:
        "0.69.0/bio/vep/plugins"


rule annotate_variants:
    input:
        calls="results/{date}/calls/ref~main/{sample}.bcf",
        plugins="resources/vep/plugins",
        fasta="resources/genomes/main.fasta",
        fai="resources/genomes/main.fasta.fai",
        gff="resources/annotation.gff.gz",
        csi="resources/annotation.gff.gz.tbi",
        problematic="resources/problematic-sites.vcf.gz",
        problematic_tbi="resources/problematic-sites.vcf.gz.tbi",
    output:
        calls="results/{date}/annotated-calls/ref~main/{sample}.bcf",
        stats="results/{date}/annotated-calls/ref~main/{sample}.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool"],
        extra=get_vep_args,
    log:
        "logs/{date}/vep/{sample}.log",
    wrapper:
        "0.72.0/bio/vep/annotate"
