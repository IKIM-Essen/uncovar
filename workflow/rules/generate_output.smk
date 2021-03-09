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
            category="3. Sequencing Details",
            subcategory="2. Read Coverage",
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


rule snakemake_report:
    input:
        "results/{date}/plots/all.strains.pangolin.svg",
        lambda wildcards: expand(
            "results/{{date}}/vcf-report/{target}.{filter}",
            target=get_samples_for_date(wildcards.date) + ["all"],
            filter=config["variant-calling"]["filters"],
        ),
        "results/{date}/plots/coverage.svg",
        "results/{date}/qc/laboratory/multiqc.html",
        "results/rki/{date}_uk-essen_rki.csv",
        "results/rki/{date}_uk-essen_rki.fasta",
        lambda wildcards: expand(
            "results/{{date}}/polished-contigs/{sample}.fasta",
            sample=get_samples_for_date(wildcards.date),
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
        "snakemake --nolock --report-stylesheet resources/custom-stylesheet.css {input} --report {output} {params.for_testing}"
