rule align_sanger:
    input:
        target="resources/genomes/main.fasta",
        query=get_sanger_files_for_sample,
    output:
        "results/{date}/sanger-aligned/ref~main/{sample}.sam",
    log:
        "results/{date}/sanger-aligned/ref~main/{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -a {input.target} {input.query} -o {output} 2> {log}"


rule sort_bam_sanger:
    input:
        "results/{date}/sanger-aligned/ref~{reference}/{sample}.sam",
    output:
        "results/{date}/sanger-aligned/ref~{reference}/{sample}_sorted.bam",
    log:
        "results/{date}/sam-to-bam/ref~{reference}/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -o {output} -O BAM {input}"


rule freebayes_sanger:
    input:
        ref="resources/genomes/main.fasta",
        ref_idx="resources/genomes/main.fasta.fai",
        # you can have a list of samples here
        samples="results/{date}/sanger-aligned/ref~main/{sample}_sorted.bam",
        index="results/{date}/sanger-aligned/ref~main/{sample}_sorted.bam.bai",
    output:
        "results/{date}/sanger-var-calls/ref~main/{sample}.bcf",
    params:
        # genotyping is performed by varlociraptor, hence we deactivate it in freebayes by
        # always setting --pooled-continuous
        extra=("--min-alternate-count 1"),
    log:
        "logs/{date}/sanger-aligned/freebayes/ref~main/{sample}.log",
    wrapper:
        "0.68.0/bio/freebayes"


rule annotate_variants_sanger:
    input:
        calls="results/{date}/sanger-var-calls/ref~main/{sample}.bcf",
        fasta="resources/genomes/main.fasta",
        fai="resources/genomes/main.fasta.fai",
        gff="resources/annotation.gff.gz",
        csi="resources/annotation.gff.gz.tbi",
        plugins="resources/vep/plugins",
        problematic="resources/problematic-sites.vcf.gz",
        problematic_tbi="resources/problematic-sites.vcf.gz.tbi",
    output:
        calls="results/{date}/sanger-var-calls_annotated/ref~main/{sample}.bcf",
        stats="results/{date}/sanger-var-calls_annotated/ref~main/{sample}.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=[],
        extra=get_vep_args,
    log:
        "logs/{date}/vep_sanger/main/{sample}.log",
    wrapper:
        "0.72.0/bio/vep/annotate"


rule compare_sanger2:
    input:
        ngs_calls="results/{date}/annotated-calls/ref~main/{sample}.bcf",
        sanger_calls="results/{date}/sanger-var-calls_annotated/ref~main/{sample}.bcf",
        sanger_aln="results/{date}/sanger-aligned/ref~main/{sample}.sam",
    output:
        sanger_vs_ngs="results/{date}/sanger-vs-genome/{sample}.txt",
    log:
        "logs/{date}/sanger-vs-genome/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/sanger-comparison2.py"


rule get_comps:
    input:
        expand(
            "results/{date}/sanger-vs-genome/{sample}.txt",
            date="2022-01-25",
            sample=get_samples(),
        ),


# rule compare_sanger:
#     input:
#         sanger=lambda wildcards: get_sanger_files(
#             wildcards,
#             "expand-sample",
#             "results/{{date}}/sanger-var-calls/ref~main/annotated_{region}~{sample}.bcf",
#         ),
#         ngs_genome="results/{date}/annotated-calls/ref~main/{sample}.bcf",
#         sanger_vs_ngs_genome=lambda wildcards: get_sanger_files(
#             wildcards,
#             "expand-sample",
#             "results/{{date}}/sanger-aligned/ref~{sample}/{region}~{sample}.csv",
#         ),
#         coverage="results/{date}/qc/samtools_depth/{sample}.txt",
#     output:
#         sanger_ngs_diff="results/{date}/sanger-vs-genome/{sample}.txt",
#         sanger_vars="results/{date}/sanger-vs-genome/vars/sanger_{sample}.csv",
#         ngs_vars="results/{date}/sanger-vs-genome/vars/ngs_{sample}.csv",
#         sanger_vs_genome="results/{date}/sanger-vs-genome/{sample}.csv",
#         sanger_ngs_diff_readable="results/{date}/sanger-vs-genome/{sample}_readable.txt",
#     log:
#         "logs/{date}/sanger-vs-genome/{sample}.log",
#     params:
#         voc=config.get("voc"),
#         regions=lambda wildcards: get_sanger_files(wildcards, "regions"),
#     script:
#         "../workflow/scripts/sanger-comp.py"


rule aggregate_sanger:
    input:
        sanger=lambda wildcards: expand(
            "results/{date}/sanger-vs-genome/vars/sanger_{sample}.csv",
            sample=get_sanger_files(wildcards, "samples"),
            date=get_latest_run_date(),
        ),
        ngs=lambda wildcards: expand(
            "results/{date}/sanger-vs-genome/vars/ngs_{sample}.csv",
            sample=get_sanger_files(wildcards, "samples"),
            date=get_latest_run_date(),
        ),
    output:
        sanger="results/{date}/sanger-vs-genome/vars/sanger_all.csv",
        ngs="results/{date}/sanger-vs-genome/vars/ngs_all.csv",
    log:
        "logs/{date}/sanger-vs-genome/all.log",
    params:
        voc=config.get("voc"),
    shell:
        "cat {input.sanger} > {output.sanger} && "
        "cat {input.ngs} > {output.ngs}"


rule count_sanger:
    input:
        lambda wildcards: expand(
            "results/{{date}}/sanger-vs-genome/{sample}.txt",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        "results/{date}/percent_sanger/sanger_all.txt",
    log:
        "logs/{date}/percent_sanger/all.log",
    script:
        "../workflow/scripts/sanger-vars.py"


rule plot_ngs_coverage_for_sanger:
    input:
        variants=expand_samples_for_date(
            "results/{{date}}/sanger-vs-genome/vars/sanger_{sample}.csv"
        ),
    output:
        table="results/{date}/plot-ngs-coverage/all.csv",
        plot="results/{date}/plot-ngs-coverage/all.svg",
    log:
        "logs/{date}/plot-ngs-coverage/all.log",
    conda:
        "../workflow/envs/python.yaml"
    script:
        "../workflow/scripts/plot-ngs-coverage.py"
