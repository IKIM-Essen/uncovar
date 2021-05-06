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


rule plot_coverage_main_sequence:
    input:
        lambda wildcards: expand(
            "results/{{date}}/qc/samtools_depth/{sample}.txt",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        report(
            "results/{date}/plots/coverage-reference-genome.svg",
            caption="../report/all-main-coverage.rst",
            category="3. Sequencing Details",
            subcategory="2. Read Coverage of Reference Genome",
        ),
    log:
        "logs/{date}/plot-coverage-main-seq.log",
    params:
        min_coverage=config["RKI-quality-criteria"]["min-depth-with-PCR-duplicates"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-all-coverage.py"


rule plot_coverage_final_sequence:
    input:
        lambda wildcards: expand(
            "results/{{date}}/tables/coverage/{sample}.txt",
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        report(
            "results/{date}/plots/coverage-assembled-genome.svg",
            caption="../report/all-final-coverage.rst",
            category="3. Sequencing Details",
            subcategory="3. Read Coverage of Reconstructed Genome",
        ),
    log:
        "logs/{date}/plot-coverage-final-seq.log",
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


rule rki_report:
    input:
        filtered_samples="results/{date}/rki-filter/{date}.txt",
        contigs=lambda wildcards: expand(
            "results/{date}/contigs-masked/{sample}.fasta",
            zip,
            date=[wildcards.date] * len(get_samples_for_date(wildcards.date)),
            sample=get_samples_for_date(wildcards.date, filtered=True),
        ),
    output:
        fasta=report(
            "results/rki/{date}_uk-essen_rki.fasta",
            category="4. RKI Submission",
            caption="../report/rki-submission-fasta.rst",
        ),
        table=report(
            "results/rki/{date}_uk-essen_rki.csv",
            category="4. RKI Submission",
            caption="../report/rki-submission-csv.rst",
        ),
    params:
        min_length=config["rki-output"]["minimum-length"],
    log:
        "logs/{date}/rki-output/{date}.log",
    script:
        "../scripts/generate-rki-output.py"


rule virologist_report:
    input:
        reads_unfiltered=lambda wildcards: expand(
            "results/{{date}}/trimmed/{sample}.fastp.json",
            sample=get_samples_for_date(wildcards.date),
        ),
        reads_used_for_assembly=lambda wildcards: expand(
            "results/{{date}}/tables/read_counts/{sample}.txt",
            sample=get_samples_for_date(wildcards.date),
        ),
        initial_contigs=lambda wildcards: get_expanded_contigs(wildcards),
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
    log:
        "logs/{date}/viro_report.log",
    params:
        voc=config.get("voc"),
        samples=lambda wildcards: get_samples_for_date(wildcards.date),
    conda:
        "../envs/pysam.yaml"
    threads: 1
    script:
        "../scripts/generate_virologist_output.py"


rule qc_html_report:
    input:
        "results/{date}/virologist/qc_report.csv",
    output:
        report(
            directory("results/{date}/qc_data/"),
            htmlindex="index.html",
            caption="../report/qc-report.rst",
            category="1. Overview",
            subcategory="1. QC Report",
        ),
    conda:
        "../envs/rbt.yaml"
    params:
        formatter=get_resource("report-table-formatter.js"),
        pin_until="sample",
    log:
        "logs/{date}/qc_report_html.log",
    shell:
        "rbt csv-report {input} --formatter {params.formatter} --pin-until {params.pin_until} {output} > {log} 2>&1"


rule snakemake_reports:
    input:
        "results/{date}/plots/coverage-reference-genome.svg",
        "results/{date}/plots/coverage-assembled-genome.svg",
        lambda wildcards: expand(
            "results/{{date}}/polished-contigs/{sample}.fasta",
            sample=get_samples_for_date(wildcards.date),
        ),
        # lambda wildcards: expand(
        #     "results/{{date}}/plots/strain-calls/{sample}.strains.kallisto.svg",
        #     sample=get_samples_for_date(wildcards.date),
        # ) if config["strain-calling"]["use-kallisto"] else "",
        "results/{date}/qc_data",
        # expand(
        #     "results/{{date}}/plots/all.{mode}-strain.strains.kallisto.svg",
        #     mode=["major", "any"],
        # ) if config["strain-calling"]["use-kallisto"] else "",
        # lambda wildcards: expand(
        #     "results/{{date}}/plots/strain-calls/{sample}.strains.pangolin.svg",
        #     sample=get_samples_for_date(wildcards.date),
        # ),
        "results/{date}/plots/all.strains.pangolin.svg",
        lambda wildcards: expand(
            "results/{{date}}/vcf-report/{target}.{filter}",
            target=get_samples_for_date(wildcards.date) + ["all"],
            filter=config["variant-calling"]["filters"],
        ),
        "results/{date}/qc/laboratory/multiqc.html",
        "results/rki/{date}_uk-essen_rki.csv",
        "results/rki/{date}_uk-essen_rki.fasta",
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
        "logs/snakemake_reports/{date}.log",
    shell:
        "snakemake --nolock --report-stylesheet resources/custom-stylesheet.css {input} --report {output} {params.for_testing}"
