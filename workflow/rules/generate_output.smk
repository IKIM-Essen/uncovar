rule masking:
    input:
        bamfile="results/{date}/mapped/ref~polished-{sample}/{sample}.bam",
        sequence="results/{date}/polished-contigs/{sample}.fasta",
    output:
        masked_sequence="results/{date}/contigs-masked/{sample}.fasta",
        coverage="results/{date}/tables/coverage/{sample}.txt",
    params:
        min_coverage=config["RKI-quality-criteria"]["min-depth-with-PCR-duplicates"],
        min_allele=config["RKI-quality-criteria"]["min-allele"],
    log:
        "logs/{date}/masking/{sample}.logs",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/mask-contigs.py"


rule plot_coverage:
    input:
        lambda wildcards: expand(
            "results/{{date}}/tables/coverage/{sample}.txt",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        report(
            "results/{date}/plots/coverage.svg",
            caption="../report/all-coverage.rst",
            category="Coverage",
        ),
    log:
        "logs/{date}/plot-coverage.log",
    params:
        min_coverage=config["RKI-quality-criteria"]["min-depth-with-PCR-duplicates"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-all-coverage.py"


checkpoint rki_filter:
    input:
        quast=lambda wildcards: expand(
            "results/{date}/quast/masked/{sample}/report.tsv",
            zip,
            date=[wildcards.date] * len(get_samples_for_date(wildcards.date)),
            sample=get_samples_for_date(wildcards.date),
        ),
        contigs=lambda wildcards: expand(
            "results/{date}/contigs-masked/{sample}.fasta",
            zip,
            date=[wildcards.date] * len(get_samples_for_date(wildcards.date)),
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        temp("results/{date}/rki-filter/{date}.txt"),
    params:
        min_identity=config["RKI-quality-criteria"]["min-identity"],
        max_n=config["RKI-quality-criteria"]["max-n"],
    log:
        "logs/{date}/rki-filter.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/rki-filter.py"


rule generate_rki:
    input:
        filtered_samples="results/{date}/rki-filter/{date}.txt",
        contigs=lambda wildcards: expand(
            "results/{date}/contigs-masked/{sample}.fasta",
            zip,
            date=[wildcards.date] * len(get_samples_for_date(wildcards.date)),
            sample=get_samples_for_date(wildcards.date, filtered=True),
        ),
    output:
        fasta="results/rki/{date}_uk-essen_rki.fasta",
        table="results/rki/{date}_uk-essen_rki.csv",
    params:
        min_length=config["rki-output"]["minimum-length"],
    log:
        "logs/{date}/rki-output/{date}.log",
    script:
        "../scripts/generate-rki-output.py"


rule generate_virologist_output:
    input:
        # reads_unfiltered=lambda wildcards: [pep.sample_table.loc[sample][["fq1", "fq2"]] for sample in get_samples_for_date(wildcards.date)],
        reads_unfiltered=lambda wildcards: expand(
            "results/{{date}}/trimmed/{sample}.fastp.json",
            sample=get_samples_for_date(wildcards.date),
        ),
        reads_filtered=lambda wildcards: expand(
            "results/{{date}}/assembly/{sample}/log",
            sample=get_samples_for_date(wildcards.date),
        ),
        initial_contigs=lambda wildcards: expand(
            "results/{{date}}/assembly/{sample}/{sample}.contigs.fa",
            sample=get_samples_for_date(wildcards.date),
        ),
        polished_contigs=lambda wildcards: expand(
            "results/{{date}}/polished-contigs/{sample}.fasta",
            sample=get_samples_for_date(wildcards.date),
        ),
        kraken=lambda wildcards: expand(
            "results/{{date}}/species-diversity/{sample}/{sample}.uncleaned.kreport2",
            sample=get_samples_for_date(wildcards.date),
        ),
        pangolin=lambda wildcards: expand(
            "results/{{date}}/tables/strain-calls/{sample}.strains.pangolin.csv",
            sample=get_samples_for_date(wildcards.date),
        ),
        bcf=lambda wildcards: expand(
            "results/{{date}}/filtered-calls/ref~main/{sample}.subclonal.high+moderate-impact.bcf",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        all_data="results/{date}/virologist/report.csv",
        qc_data="results/{date}/virologist/qc_report.csv",
        var_data="results/{date}/virologist/var_report.csv",
    log:
        "logs/{date}/viro_report.log",
    params:
        voc=config.get("voc"),
    conda:
        "../envs/pysam.yaml"
    threads: 1
    script:
        "../scripts/generate_virologist_output.py"


rule snakemake_html_report_qc:
    input:
        "results/{date}/virologist/qc_report.csv",
    output:
        report(
            directory("results/{date}/qc_data/"),
            htmlindex="index.html",
            caption="../report/qc-report.rst",
            category="QC report overview",
        ),
    conda:
        "../envs/rbt.yaml"
    log:
        "logs/{date}/qc_report_html.log",
    shell:
        "rbt csv-report {input} {output} > {log} 2>&1"


rule snakemake_html_report_variants:
    input:
        "results/{date}/virologist/var_report.csv",
    output:
        report(
            directory("results/{date}/var_data/"),
            htmlindex="index.html",
            caption="../report/var-report.rst",
            category="Variant report overview",
        ),
    conda:
        "../envs/rbt.yaml"
    log:
        "logs/{date}/var_report_html.log",
    shell:
        "rbt csv-report {input} {output} > {log} 2>&1"


rule snakemake_reports:
    input:
        "results/{date}/plots/coverage.svg",
        lambda wildcards: expand(
            "results/{{date}}/polished-contigs/{sample}.fasta",
            sample=get_samples_for_date(wildcards.date),
        ),
        lambda wildcards: expand(
            "results/{{date}}/plots/strain-calls/{sample}.strains.kallisto.svg",
            sample=get_samples_for_date(wildcards.date),
        ),
        "results/{date}/qc_data",
        "results/{date}/var_data",
        expand(
            "results/{{date}}/plots/all.{mode}-strain.strains.kallisto.svg",
            mode=["major", "any"],
        ),
        lambda wildcards: expand(
            "results/{{date}}/plots/strain-calls/{sample}.strains.pangolin.svg",
            sample=get_samples_for_date(wildcards.date),
        ),
        "results/{date}/plots/all.strains.pangolin.svg",
        lambda wildcards: expand(
            "results/{{date}}/scenarios/{sample}.yaml",
            sample=get_samples_for_date(wildcards.date),
        ),
        lambda wildcards: expand(
            "results/{{date}}/vcf-report/{target}.{filter}",
            target=get_samples_for_date(wildcards.date) + ["all"],
            filter=config["variant-calling"]["filters"],
        ),
        expand(
            "results/{{date}}/ucsc-vcfs/all.{{date}}.{filter}.vcf",
            filter=config["variant-calling"]["filters"],
        ),
    output:
        "results/reports/{date}.zip",
    params:
        for_testing=(
            "--snakefile ../workflow/Snakefile"
            if config.get("benchmark-genomes", [])
            else ""
        ),
    log:
        "../logs/snakemake_reports/{date}.log",
    shell:
        "snakemake --nolock {input} --report {output} {params.for_testing}"
