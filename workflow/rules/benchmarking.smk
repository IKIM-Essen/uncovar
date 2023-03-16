# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule simulate_strain_reads:
    input:
        get_genome_fasta,
    output:
        left=temp("resources/benchmarking/{accession}/reads.1.fastq.gz"),
        right=temp("resources/benchmarking/{accession}/reads.2.fastq.gz"),
    params:
        no_reads=lambda wildcards: no_reads(wildcards),
        length_reads=lambda wildcards: length_read(wildcards),
    log:
        "logs/mason/benchmarking/{accession}.log",
    conda:
        "../envs/mason.yaml"
    threads: 4
    shell:  # median reads in data: 584903
        "mason_simulator -ir {input} -n {params.no_reads} --illumina-read-length {params.length_reads} --num-threads {threads} -o {output.left} -or {output.right} --fragment-mean-size 400 2> {log}"


rule mix_strain_reads:
    input:
        left=expand(
            "resources/benchmarking/{mix}/reads.1.fastq.gz",
            mix=[
                "{MIXTURE_PART_INDICATOR}{{strain_{no}}}".format(
                    MIXTURE_PART_INDICATOR=MIXTURE_PART_INDICATOR, no=i
                )
                for i in range(config["mixtures"]["no_strains"])
            ],
        ),
        right=expand(
            "resources/benchmarking/{mix}/reads.2.fastq.gz",
            mix=[
                "{MIXTURE_PART_INDICATOR}{{strain_{no}}}".format(
                    MIXTURE_PART_INDICATOR=MIXTURE_PART_INDICATOR, no=i
                )
                for i in range(config["mixtures"]["no_strains"])
            ],
        ),
    output:
        left=temp(
            expand(
                "resources/mixtures/{mix}/reads.1.fastq.gz",
                mix="".join(
                    [
                        "{MIXTURE_PART_INDICATOR}{{strain_{no}}}".format(
                            MIXTURE_PART_INDICATOR=MIXTURE_PART_INDICATOR,
                            no=i,
                        )
                        for i in range(config["mixtures"]["no_strains"])
                    ]
                ),
            )
        ),
        right=temp(
            expand(
                "resources/mixtures/{mix}/reads.2.fastq.gz",
                mix="".join(
                    [
                        "{MIXTURE_PART_INDICATOR}{{strain_{no}}}".format(
                            MIXTURE_PART_INDICATOR=MIXTURE_PART_INDICATOR,
                            no=i,
                        )
                        for i in range(config["mixtures"]["no_strains"])
                    ]
                ),
            )
        ),
    log:
        "logs/mix_strain_reads/{}".format(
            "".join(
                [
                    "{MIXTURE_PART_INDICATOR}{{strain_{no}}}".format(
                        MIXTURE_PART_INDICATOR=MIXTURE_PART_INDICATOR, no=i
                    )
                    for i in range(config["mixtures"]["no_strains"])
                ]
            )
        ),
    shell:
        "(zcat {input.left} > {output.left} &&"
        "zcat {input.right} > {output.right}) 2>{log}"


rule test_benchmark_results:
    input:
        get_benchmark_results,
    output:
        "results/benchmarking/strain-calling.csv",
    params:
        true_accessions=get_strain_accessions,
    log:
        "logs/test-benchmark-results.log",
    resources:
        notebooks=1,
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/test-benchmark-results.py.ipynb"


rule test_assembly_results:
    input:
        "resources/genomes/{accession}.fasta",
        get_assembly_result,
    output:
        "results/benchmarking/assembly/{assembly_type}/{accession}.bam",
    log:
        "logs/test-assembly-results/{assembly_type}/{accession}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 --MD --eqx -ax asm5 {input} -o {output} 2> {log}"


rule summarize_assembly_results:
    input:
        bams=get_assembly_comparisons(bams=True),
        refs=get_assembly_comparisons(bams=False),
    output:
        "results/benchmarking/assembly/{assembly_type}.csv",
    log:
        "logs/summarize-assembly-results/{assembly_type}/assembly-results.log",
    resources:
        notebooks=1,
    conda:
        "../envs/pysam.yaml"
    notebook:
        "../notebooks/assembly-benchmark-results.py.ipynb"


rule test_non_cov2:
    input:
        pangolin=get_non_cov2_calls(from_caller="pangolin"),
        kallisto=get_non_cov2_calls(from_caller="kallisto"),
    output:
        "results/benchmarking/non-sars-cov-2.csv",
    params:
        accessions=get_non_cov2_accessions(),
    log:
        "logs/benchmarking/summarize_non_cov2.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize-non-cov2.py"


rule report_non_cov2:
    input:
        summary="results/benchmarking/non-sars-cov-2.csv",
        call_plots_kallisto=expand(
            "results/benchmarking/plots/strain-calls/non-cov2-{accession}.strains.kallisto.svg",
            accession=get_non_cov2_accessions(),
        ),
        call_plots_pangolin=expand(
            "results/benchmarking/plots/strain-calls/non-cov2-{accession}.polished.strains.pangolin.svg",
            accession=get_non_cov2_accessions(),
        ),
    output:
        report(
            directory("results/benchmarking/html"),
            htmlindex="index.html",
            category="Test results",
        ),
    log:
        "logs/report_non_cov2.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report -s '\t' {input.summary} {output}"


checkpoint generate_mixtures:
    input:
        "results/benchmarking/tables/strain-genomes.txt",
    output:
        "results/benchmarking/tables/mixtures.txt",
    params:
        mixtures=generate_mixtures,
    log:
        "logs/generate_mixtures.log",
    script:
        "../scripts/generate-mixtures.py"


rule evaluate_strain_call_error:
    input:
        get_mixture_results,
    output:
        "results/benchmarking/tables/{caller}-strain-call-error.csv",
    params:
        max_reads=config["mixtures"]["max_reads"],
        prefix=MIXTURE_PREFIX,
        separator=MIXTURE_PART_INDICATOR,
        percentage=MIXTURE_PERCENTAGE_INDICATOR,
    log:
        "logs/evaluate-{caller}-strain-call-error.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/evaluate-strain-call-error.py"


rule plot_strain_call_error:
    input:
        "results/benchmarking/tables/{caller}-strain-call-error.csv",
    output:
        report(
            "results/benchmarking/plots/{caller}-strain-call-error-heatmap.svg",
            category="Figure 4: Strain Call Error",
            caption="../report/publication-strain-call-error.rst",
        ),
        "results/benchmarking/plots/{caller}-strain-call-error-false-predictions.svg",
        "results/benchmarking/plots/{caller}-strain-call-error-content-false-predictions.svg",
    log:
        "logs/plot-{caller}-strain-call-error.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-caller-error.py"


rule assembly_comparison_trinity:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        temp("results/{date}/assembly/{sample}/trinity-pe/{sample}.contigs.fasta"),
    log:
        "logs/{date}/trinity/{sample}.log",
    params:
        extra="",
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/trinity.yaml"
    shell:
        "(Trinity --left {input.fastq1} --max_memory 16G --right {input.fastq2} --CPU {threads} --seqType fq --output {params.outdir} && "
        "mv {params.outdir}.Trinity.fasta {output} ) > {log} 2>&1"


rule assembly_comparison_velvet:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        temp("results/{date}/assembly/{sample}/velvet-pe/{sample}.contigs.fasta"),
    log:
        "logs/{date}/velvet/{sample}.log",
    params:
        extra="",
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/velvet.yaml"
    shell:
        """
        velveth {params.outdir} 21 -fastq.gz -shortPaired {input.fastq1} {input.fastq2} > {log} 2>&1
        velvetg {params.outdir} -ins_length 150 -exp_cov 10 >> {log} 2>&1
        mv {params.outdir}/contigs.fa {output} >> {log} 2>&1
        """


rule order_contigs_assembly_comparison:
    input:
        contigs=(
            "results/{date}/assembly/{sample}/{assembler}-pe/{sample}.contigs.fasta"
        ),
        reference="resources/genomes/main.fasta",
    output:
        temp(
            "results/{date}/assembly/{sample}/{assembler}/{sample}.ordered.contigs.fasta"
        ),
    log:
        "logs/{date}/ragtag/{assembler}/{sample}.log",
    params:
        outdir=get_output_dir,
    conda:
        "../envs/ragtag.yaml"
    shadow:
        "minimal"
    shell:
        "(cd {params.outdir} &&"
        " ragtag.py scaffold -C ../../../../../{input.reference} ../../../../../{input.contigs} &&"
        " cd ../../../../../ && mv {params.outdir}/ragtag_output/ragtag.scaffold.fasta {output})"
        " > {log} 2>&1"


use rule filter_chr0 as filter_chr0_assembly_comparison with:
    input:
        "results/{date}/assembly/{sample}/{assembler}/{sample}.ordered.contigs.fasta",
    output:
        "results/{date}/assembly/{sample}/{assembler}/{sample}.contigs.ordered.filtered.fasta",
    log:
        "logs/{date}/ragtag/{assembler}/{sample}_cleaned.log",


use rule align_contigs as align_contigs_assembly_comparison with:
    input:
        target="resources/genomes/main.fasta",
        query="results/{date}/assembly/{sample}/{assembler}-pe/{sample}.contigs.fasta",
    output:
        "results/{date}/assembly/{sample}/{assembler}/main_{sample}.bam",
    log:
        "results/{date}/assembly/{sample}/{assembler}/main_{sample}.log",


use rule quast as quast_assembly_comparison with:
    input:
        fasta="results/{date}/assembly/{sample}/{assembler}-pe/{sample}.contigs.fasta",
        bam="results/{date}/assembly/{sample}/{assembler}/main_{sample}.bam",
        reference="resources/genomes/main.fasta",
    output:
        "results/{date}/assembly/{sample}/{assembler}/quast/report.tsv",
        "results/{date}/assembly/{sample}/{assembler}/quast/transposed_report.tsv",
    log:
        "logs/{date}/assembly/quast/{assembler}/{sample}.log",


rule plot_assemblies:
    input:
        initial=get_samples_for_assembler_comparison(
            "results/{zip1}/assembly/{zip2}/{{exp}}-pe/{zip2}.contigs.fasta"
        ),
        final=get_samples_for_assembler_comparison(
            "results/{zip1}/assembly/{zip2}/{{exp}}/{zip2}.contigs.ordered.filtered.fasta",
        ),
        quast=get_samples_for_assembler_comparison(
            "results/{zip1}/assembly/{zip2}/{{exp}}/quast/transposed_report.tsv",
        ),
    output:
        report(
            "results/benchmarking/plots/assembler-comparison.svg",
            category="Figure 2: Assembler Comparison",
            caption="../report/publication-assembler-comparison.rst",
        ),
        "results/benchmarking/plots/assembler-comparison.csv",
        report(
            "results/benchmarking/plots/assembler-comparison_genome_fraction.svg",
            category="Supplementary Figure 2: Genome Fraction",
            caption="../report/publication-genome-fraction.rst",
        ),
    log:
        "logs/benchmarking/all_assemblies_plot.log",
    params:
        samples=get_samples(),
        assembler=config["assemblers_for_comparison"],
        amplicon_state=lambda wildcards: get_list_of_amplicon_states(wildcards),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-assembly-comparison.py"


rule get_read_length_statistics:
    input:
        expand(
            "results/{date}/tables/read_pair_counts/{sample}.txt",
            zip,
            date=get_dates(),
            sample=get_samples(),
        ),
    output:
        "results/benchmarking/tables/read_statistics.txt",
    log:
        "logs/get_read_statistics.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get-read-statistics.py"


rule plot_dependency_of_pangolin_call:
    input:
        get_mixture_results,
    output:
        report(
            "results/benchmarking/plots/{caller}-call-dependency.svg",
            category="Figure 3: Lineage Call Dependency",
            caption="../report/publication-lineage-call-dependency.rst",
        ),
    log:
        "logs/plot_dependency_of_{caller}_call.log",
    params:
        prefix=MIXTURE_PREFIX,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-dependency-of-pangolin-call.py"


rule plot_pangolin_conflict:
    input:
        get_mixture_results,
    output:
        "results/benchmarking/plots/{caller}_statistics.svg",
        "results/benchmarking/tables/{caller}_statistics.csv",
    log:
        "logs/plot_pangolin_conflict_{caller}.log",
    params:
        separator=MIXTURE_PART_INDICATOR,
        percentage=MIXTURE_PERCENTAGE_INDICATOR,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-pangolin-conflict.py"


rule collect_lineage_calls_of_various_stages:
    input:
        kallisto=expand(
            "results/benchmarking/tables/strain-calls/{prefix}{{lineage}}{number_indi}{{number}}{length_indi}{{length}}{state_indi}reads.strains.kallisto.tsv",
            prefix=READ_TEST_PREFIX,
            number_indi=READ_NUMBER_INDICATOR,
            length_indi=READ_LENGTH_INDICATOR,
            state_indi=READ_STATE_INDICATOR,
        ),
        pangolin=expand(
            "results/benchmarking/tables/strain-calls/{prefix}{{lineage}}{number_indi}{{number}}{length_indi}{{length}}{state_indi}{state}.polished.strains.pangolin.csv",
            prefix=READ_TEST_PREFIX,
            number_indi=READ_NUMBER_INDICATOR,
            length_indi=READ_LENGTH_INDICATOR,
            state_indi=READ_STATE_INDICATOR,
            state=["contig", "scaffold", "polished_scaffold", "pseudo"],
        ),
    output:
        "results/benchmarking/tables/collected_lineage_calls_on_{lineage}_{number}_{length}.tsv",
    params:
        states=["contig", "scaffold", "polished_scaffold", "pseudo"],
    log:
        "logs/collect_lineage_calls/{lineage}_{number}_{length}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/collect_lineage_calls.py"


rule get_largest_contig:
    input:
        "results/{date}/assembly/megahit/{sample}/{sample}.contigs.fasta",
    output:
        "results/{date}/tables/largest_contig/{sample}.fasta",
    log:
        "logs/{date}/get_largest_contig/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_largest_contig.py"


checkpoint select_random_lineages:
    input:
        "results/{date}/tables/strain-genomes.txt",
    output:
        "results/{date}/tables/selected-strain-genomes-reads.txt",
    params:
        number_of_samples=config["read_lineage_call"]["number_of_samples"],
    log:
        "logs/{date}/select_random_lineages.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/select_random_lineages.py"


rule aggregate_read_calls:
    input:
        get_read_calls,
    output:
        "results/benchmarking/tables/aggregated_read_calls.tsv",
    log:
        "logs/aggregate_read_calls.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/aggregate_read_calls.py"


rule plot_read_call:
    input:
        "results/benchmarking/tables/aggregated_read_calls.tsv",
    output:
        "results/benchmarking/plots/aggregated_read_calls.svg",
    log:
        "logs/plot_read_call.log",
    resources:
        notebooks=1,
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/plot-read-call.py.ipynb"


rule get_publication_plots:
    input:
        expand(
            [
                "results/benchmarking/plots/{caller}-strain-call-error-heatmap.svg",
                "results/benchmarking/plots/{caller}-call-dependency.svg",
            ],
            caller=["kallisto", "pangolin"],
        ),
        "results/benchmarking/plots/assembler-comparison.svg",
        "results/benchmarking/plots/assembler-comparison_genome_fraction.svg",


rule indentify_test_case_variants:
    input:
        illumina_bcf=get_test_cases_variant_calls(ILLUMINA),
        illumina_csi=get_test_cases_variant_calls(ILLUMINA, ".csi"),
        ont_bcf=get_test_cases_variant_calls(ONT),
        ont_csi=get_test_cases_variant_calls(ONT, ".csi"),
    output:
        "results/benchmarking/tables/test-cases/{test_case}/found-variants.tsv",
    log:
        "logs/test_cases/{test_case}.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/benchmarking/compare-vcf.py"


rule aggregate_test_case_variants:
    input:
        lambda wildcards: expand(
            "results/benchmarking/tables/test-cases/{test_case}/found-variants.tsv",
            test_case=get_all_test_cases_names(wildcards),
        ),
    output:
        "results/benchmarking/tables/test-cases/found-variants.tsv",
    log:
        "logs/aggregate_test_case_variants.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/benchmarking/aggregate-test-case-variants.py"


rule filter_test_case_variants:
    input:
        "results/benchmarking/tables/test-cases/found-variants.tsv",
    output:
        different_probs="results/benchmarking/tables/test-cases/different-prob-calls.tsv",
        illumina_only="results/benchmarking/tables/test-cases/illumina-only-calls.tsv",
        ont_only="results/benchmarking/tables/test-cases/ont-only-calls.tsv",
    log:
        "logs/filter_test_case_variants.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/benchmarking/filter-test-case-variants.py"


checkpoint get_test_case_variant_paths:
    input:
        "results/benchmarking/tables/test-cases/different-prob-calls.tsv",
    output:
        overview="results/testcases/different-probs-overview.tsv",
        paths="results/benchmarking/tables/test-cases/aggregated-variants.tsv",
    params:
        illumina=ILLUMINA,
        ont=ONT,
        illumina_varrange=ILLUMINA_VARRANGE,
        ont_varrange=ONT_VARRANGE,
        sample_table=pep.sample_table.to_dict(),
    log:
        "logs/get_test_case_variant_paths.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/benchmarking/get-test-case-variant-paths.py"


checkpoint check_presence_of_test_case_variant_in_call:
    input:
        bcfs=get_aggregated_test_case_variants("bcf"),
        csi=get_aggregated_test_case_variants("csi"),
    output:
        "results/benchmarking/tables/test-cases/presence-test-case-varaints.tsv",
    params:
        variants=get_aggregated_test_case_variants("variants"),
        poses=get_aggregated_test_case_variants("poses"),
        test_cases=get_aggregated_test_case_variants("test-case-paths"),
    log:
        "logs/check_presence_of_test_case_variant.log",
    conda:
        "../envs/pysam.yaml"
    script:
        "../scripts/benchmarking/check-presence-of-test-case-variant-in-call.py"


rule varlociraptor_test_case:
    input:
        obs="results/{date}/observations/ref~{reference}/{sample}.{varrange}.bcf",
        scenario="results/{date}/scenarios/{sample}.yaml",
    output:
        bcf=temp(
            "results/{date}/call-test-cases/ref~{reference}/{sample}.{varrange}.chrom~{chrom}.pos~{pos}.bcf"
        ),
        testcase=directory(
            "results/testcases/{sample}@{chrom}:{pos}-from:{date}-ref:{reference}-varrange:{varrange}."
        ),
    params:
        biases=get_varlociraptor_bias_flags,
    log:
        "logs/{date}/varlociraptor/call/ref~{reference}/chrom~{chrom}.pos~{pos}.{sample}.{varrange}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor call variants "
        "--testcase-prefix {output.testcase} "
        "--testcase-locus {wildcards.chrom}:{wildcards.pos} "
        "{params.biases} generic --obs {wildcards.sample}={input.obs} "
        "--scenario {input.scenario} > {output.bcf} 2> {log}"
