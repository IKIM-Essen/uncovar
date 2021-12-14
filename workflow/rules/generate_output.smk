# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule masking:
    input:
        bamfile="results/{date}/mapped/ref~{reference}-{sample}/{sample}.bam",
        bai="results/{date}/mapped/ref~{reference}-{sample}/{sample}.bam.bai",
        sequence="results/{date}/contigs/{reference}/{sample}.fasta",
    output:
        masked_sequence="results/{date}/contigs/masked/{reference}/{sample}.fasta",
        coverage="results/{date}/tables/coverage/{reference}/{sample}.txt",
    params:
        min_coverage=config["quality-criteria"]["min-depth-with-PCR-duplicates"],
        min_allele=config["quality-criteria"]["min-allele"],
        is_ont=is_ont,
    log:
        "logs/{date}/masking/{reference}/{sample}.logs",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/masking.py"


rule plot_coverage_main_sequence:
    input:
        expand_samples_for_date("results/{{date}}/qc/samtools_depth/{sample}.txt"),
    output:
        report(
            "results/{date}/plots/coverage-reference-genome.svg",
            caption="../report/all-main-coverage.rst",
            category="3. Sequencing Details",
            subcategory="2. Coverage of Reference Genome",
        ),
    params:
        min_coverage=config["quality-criteria"]["min-depth-with-PCR-duplicates"],
    log:
        "logs/{date}/plot-coverage-main-seq.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-all-coverage.py"


rule plot_coverage_polished_sequence:
    input:
        expand_samples_for_date(
            "results/{{date}}/tables/coverage/polished/{sample}.txt"
        ),
    output:
        report(
            "results/{date}/plots/coverage-assembled-genome.svg",
            caption="../report/all-final-coverage.rst",
            category="3. Sequencing Details",
            subcategory="3. Coverage of Reconstructed Sequences",
        ),
    params:
        min_coverage=config["quality-criteria"]["min-depth-with-PCR-duplicates"],
    log:
        "logs/{date}/plot-coverage-final-seq.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-all-coverage.py"


checkpoint quality_filter:
    input:
        quast=get_final_assemblies_identity,
        contigs=get_final_assemblies,
    output:
        passed_filter="results/{date}/tables/quality-filter/{assembly_type}.txt",
        filter_summary="results/{date}/tables/filter_summary/{assembly_type}.tsv",
    params:
        min_identity=config["quality-criteria"]["min-identity"],
        max_n=config["quality-criteria"]["max-n"],
    log:
        "logs/{date}/quality-filter/{assembly_type}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/quality-filter.py"


rule high_quality_genomes_report:
    input:
        contigs=lambda wildcards: get_assemblies_for_submission(
            wildcards, "accepted samples"
        ),
    output:
        fasta=report(
            "results/high-quality-genomes/{date}.fasta",
            category="6. High Quality Genomes",
            caption="../report/rki-submission-fasta.rst",
        ),
        table=report(
            "results/high-quality-genomes/{date}.csv",
            category="6. High Quality Genomes",
            caption="../report/rki-submission-csv.rst",
        ),
    params:
        seq_type=lambda wildcards: get_assemblies_for_submission(
            wildcards, "accepted samples technology"
        ),
    conda:
        "../envs/pysam.yaml"
    log:
        "logs/{date}/high_quality_output.log",
    script:
        "../scripts/generate-high-quality-report.py"


rule overview_table_csv:
    input:
        reads_raw=get_raw_reads_counts,
        reads_trimmed=get_trimmed_reads_counts,
        reads_used_for_assembly=expand_samples_for_date(
            "results/{{date}}/tables/read_pair_counts/{sample}.txt",
        ),
        initial_contigs=expand_samples_for_date(
            "results/{{date}}/contigs/checked/{sample}.fasta",
        ),
        polished_contigs=expand_samples_for_date(
            "results/{{date}}/contigs/masked/polished/{sample}.fasta",
        ),
        pseudo_contigs=get_fallbacks_for_report("pseudo"),
        consensus_contigs=get_fallbacks_for_report("consensus"),
        kraken=get_kraken_output,
        pangolin=expand_samples_for_date(
            "results/{{date}}/tables/strain-calls/{sample}.polished.strains.pangolin.csv",
        ),
        bcf=expand_samples_for_date(
            "results/{{date}}/filtered-calls/ref~main/{sample}.subclonal.high+moderate-impact.bcf",
        ),
        # Added because WorkflowError: Rule parameter depends on checkpoint but checkpoint output is not defined 
        # as input file for the rule. Please add the output of the respective checkpoint to the rule inputs.
        _= expand("results/{{date}}/tables/quality-filter/{assembly_type}.txt", assembly_type=["masked-assembly", "pseudo-assembly"])
    output:
        qc_data="results/{date}/tables/overview.csv",
    params:
        assembly_used=lambda wildcards: get_assemblies_for_submission(
            wildcards, "all samples"
        ),
        voc=config.get("voc"),
        samples=lambda wildcards: get_samples_for_date(wildcards.date),
    log:
        "logs/{date}/overview-table.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/generate-overview-table.py"


rule overview_table_html:
    input:
        "results/{date}/tables/overview.csv",
    output:
        report(
            directory("results/{date}/overview/"),
            htmlindex="index.html",
            caption="../report/qc-report.rst",
            category="1. Overview",
            subcategory="1. Report",
        ),
    params:
        formatter=get_resource("report-table-formatter.js"),
        pin_until="Sample",
    log:
        "logs/{date}/qc_report_html.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report {input} --formatter {params.formatter} --pin-until {params.pin_until} {output} > {log} 2>&1"


rule filter_overview:
    input:
        de_novo="results/{date}/tables/filter_summary/masked-assembly.tsv",
        pseudo=get_if_any_pseudo_assembly(
            "results/{date}/tables/filter_summary/pseudo-assembly.tsv"
        ),
        consensus=get_if_any_consensus_assembly(
            "results/{date}/tables/filter_summary/consensus-assembly.tsv"
        ),
    output:
        "results/{date}/tables/filter-overview.csv",
    params:
        samples=lambda wildcards: get_samples_for_date(wildcards.date),
        min_identity=config["quality-criteria"]["min-identity"],
        max_n=config["quality-criteria"]["max-n"],
    log:
        "logs/{date}/filter-overview.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/generate-filter-overview.py"


rule filter_overview_html:
    input:
        "results/{date}/tables/filter-overview.csv",
    output:
        report(
            directory("results/{date}/filter-overview"),
            htmlindex="index.html",
            caption="../report/filter-overview.rst",
            category="4. Assembly",
            subcategory="0. Quality Overview",
        ),
    params:
        pin_until="Sample",
    log:
        "logs/{date}/filter_overview_html.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report {input} --pin-until {params.pin_until} {output} > {log} 2>&1"


rule plot_lineages_over_time:
    input:
        lambda wildcards: expand(
            "results/{date}/tables/strain-calls/{sample}.polished.strains.pangolin.csv",
            zip,
            date=get_dates_before_date(wildcards),
            sample=get_samples_before_date(wildcards),
        ),
    output:
        report(
            "results/{date}/plots/lineages-over-time.svg",
            caption="../report/lineages-over-time.rst",
            category="1. Overview",
            subcategory="2. Lineages Development",
        ),
        "results/{date}/tables/lineages-over-time.csv",
    params:
        dates=get_dates_before_date,
    log:
        "logs/{date}/plot_lineages_over_time.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-lineages-over-time.py"


rule pangolin_call_overview:
    input:
        get_aggregated_pangolin_calls,
    output:
        "results/{date}/tables/pangolin_calls_per_stage.csv",
    params:
        samples=lambda wildcards: get_aggregated_pangolin_calls(
            wildcards, return_list="samples"
        ),
        stages=lambda wildcards: get_aggregated_pangolin_calls(
            wildcards, return_list="stages"
        ),
    log:
        "logs/{date}/aggregate_pangolin_calls.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/aggregate-pangolin-calls-per-stage.py"


rule pangolin_call_overview_html:
    input:
        "results/{date}/tables/pangolin_calls_per_stage.csv",
    output:
        report(
            directory("results/{date}/pangolin-call-overview"),
            htmlindex="index.html",
            caption="../report/pangolin-call-overview.rst",
            category="4. Assembly",
            subcategory="0. Quality Overview",
        ),
    params:
        pin_until="Sample",
    log:
        "logs/{date}/pangolin_call_overview_html.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report {input} --pin-until {params.pin_until} {output} > {log} 2>&1"


rule snakemake_reports:
    input:
        # 1. Overview
        "results/{date}/overview/",
        "results/{date}/plots/lineages-over-time.svg",
        "results/{date}/plots/all.strains.pangolin.svg",
        "results/{date}/plots/all.major-strain.strains.kallisto.svg",
        expand_samples_for_date(
            ["results/{{date}}/plots/strain-calls/{sample}.strains.kallisto.svg"]
        ),
        # 2. Variant Call Details
        lambda wildcards: expand(
            "results/{{date}}/vcf-report/{target}.{filter}",
            target=get_samples_for_date(wildcards.date) + ["all"],
            filter=config["variant-calling"]["filters"],
        ),
        # 3. Sequencing Details
        "results/{date}/qc/laboratory/multiqc.html",
        "results/{date}/plots/coverage-reference-genome.svg",
        "results/{date}/plots/coverage-assembled-genome.svg",
        lambda wildcards: "results/{date}/plots/primer-clipping-intervals.svg"
        if any_sample_is_amplicon(wildcards)
        else [],
        # 4. Assembly
        "results/{date}/filter-overview",
        "results/{date}/pangolin-call-overview",
        expand_samples_for_date(
            [
                "results/{{date}}/contigs/polished/{sample}.fasta",
                "results/{{date}}/contigs/fallback/{sample}.fasta",
            ]
        ),
        # 5. Variant Call Files
        expand(
            ["results/{{date}}/ucsc-vcfs/all.{{date}}.{filter}.vcf"],
            filter=config["variant-calling"]["filters"],
        ),
        # 6. High Quality Genomes
        "results/high-quality-genomes/{date}.fasta",
        "results/high-quality-genomes/{date}.csv",
    output:
        "results/reports/{date}.zip",
    params:
        for_testing=get_if_testing("--snakefile ../workflow/Snakefile"),
    conda:
        "../envs/snakemake.yaml"
    log:
        "logs/snakemake_reports/{date}.log",
    shell:
        "snakemake --nolock --report-stylesheet resources/custom-stylesheet.css {input} "
        "--report {output} {params.for_testing} "
        "> {log} 2>&1"
