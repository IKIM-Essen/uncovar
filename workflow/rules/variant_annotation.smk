rule get_vep_plugins:
    output:
        directory("resources/vep/plugins"),
    params:
        release=100,
    log:
        "logs/vep-plugins.log",
    wrapper:
        "0.69.0/bio/vep/plugins"


# takes all given information and adds the genome features to the variant calling; estimates effects on proteins (indel etc)
rule annotate_variants:
    input:
        calls="results/calls/ref~main/{sample}.bcf",
        plugins="resources/vep/plugins",
        fasta="resources/genomes/main.fasta",
        fai="resources/genomes/main.fasta.fai",
        gff="resources/annotation.gff.gz",
        csi="resources/annotation.gff.gz.tbi",
        problematic="resources/problematic-sites.vcf.gz",
        problematic_tbi="resources/problematic-sites.vcf.gz.tbi",
    output:
        calls="results/annotated-calls/ref~main/{sample}.bcf",
        stats="results/annotated-calls/ref~main/{sample}.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=["LoFtool"],
        extra=lambda w, input: "--vcf_info_field ANN --hgvsg --hgvs --synonyms {synonyms} --custom {input.problematic},,vcf,exact,0,".format(
            input=input, synonyms=get_resource("synonyms.txt")
        ),
    log:
        "logs/vep/{sample}.log",
    threads: 4
    wrapper:
        "vep-nocache/bio/vep/annotate"
