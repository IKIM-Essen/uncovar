rule align_sanger:
    input:
        target="resources/genomes/main.fasta",
        query=get_aln_fastas,
    output:
        "results/{date}/sanger-aligned/ref~main/{region}~{sample}.sam",
    log:
        "results/{date}/sanger-aligned/ref~main/{region}~{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax splice {input.target} {input.query} -o {output} 2> {log}"


rule sort_sam:
    input:
        "results/{date}/sanger-aligned/ref~main/{region}~{sample}.sam",
    output:
        "results/{date}/sanger-aligned/ref~main/{region}~{sample}.bam",
    log:
        "results/{date}/sam-to-bam/ref~main/{region}~{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -o {output} -O BAM {input}"


rule freebayes_sanger:
    input:
        ref="resources/genomes/main.fasta",
        ref_idx="resources/genomes/main.fasta.fai",
        # you can have a list of samples here
        samples="results/{date}/sanger-aligned/ref~main/{region}~{sample}.bam",
        index="results/{date}/sanger-aligned/ref~main/{region}~{sample}.bam.bai",
    output:
        "results/{date}/sanger-var-calls/ref~main/{region}~{sample}.vcf",
    params:
        # genotyping is performed by varlociraptor, hence we deactivate it in freebayes by 
        # always setting --pooled-continuous
        extra=("--min-alternate-count 1"),
    log:
        "logs/{date}/sanger-aligned/freebayes/ref~main/{region}~{sample}.log",
    wrapper:
        "0.68.0/bio/freebayes"


rule annotate_variants_sanger:
    input:
        calls="results/{date}/sanger-var-calls/ref~main/{region}~{sample}.vcf",
        fasta="resources/genomes/main.fasta",
        fai="resources/genomes/main.fasta.fai",
        gff="resources/annotation.gff.gz",
        csi="resources/annotation.gff.gz.tbi",
        plugins="resources/vep/plugins",
        problematic="resources/problematic-sites.vcf.gz",
        problematic_tbi="resources/problematic-sites.vcf.gz.tbi",
    output:
        calls="results/{date}/sanger-var-calls/ref~main/annotated_{region}~{sample}.bcf",
        stats="results/{date}/sanger-var-calls/ref~main/annotated_{region}~{sample}.html",
    params:
        # Pass a list of plugins to use, see https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html
        # Plugin args can be added as well, e.g. via an entry "MyPlugin,1,FOO", see docs.
        plugins=[],
        extra=get_vep_args,
    log:
        "logs/{date}/vep_sanger/{region}~{sample}.log",
    wrapper:
        "0.72.0/bio/vep/annotate"


rule compare_sanger:
    input:
        sanger=get_sanger_files,
        genome="results/{date}/annotated-calls/ref~main/normed_{sample}.bcf",
    output:
        "results/{date}/sanger-vs-genome/{sample}.txt",
        "results/{date}/vars/sanger_{sample}.csv",
        "results/{date}/vars/ngs_{sample}.csv",
    log:
        "logs/{date}/sanger-vs-genome/{sample}.log",
    params:
        voc=config.get("voc"),
    script:
        "../scripts/sanger-comp.py"


rule aggregate_sanger:
    input:
        sanger=lambda wildcards: expand(
            "results/{date}/vars/sanger_{sample}.csv",
            sample=get_sanger_files(wildcards, "samples"),
            date=get_latest_run_date(),
        ),
        ngs=lambda wildcards: expand(
            "results/{date}/vars/ngs_{sample}.csv",
            sample=get_sanger_files(wildcards, "samples"),
            date=get_latest_run_date(),
        ),
    output:
        sanger="results/{date}/vars/sanger_all.csv",
        ngs="results/{date}/vars/ngs_all.csv",
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
        "../scripts/sanger-vars.py"


rule plot_ngs_coverage_for_sanger:
    input:
        coverage=expand_samples_for_date(
            "results/{{date}}/qc/samtools_depth/{sample}.txt"
        ),
        variants=expand_samples_for_date("results/{{date}}/vars/ngs_{sample}.csv"),
    output:
        table="results/{date}/plot-ngs-coverage/all.csv",
        # plot="results/{date}/plot-ngs-coverage/all.svg",
    log:
        "logs/{date}/plot-ngs-coverage/all.log",
    script:
        "../scripts/plot-ngs-coverage.py"
