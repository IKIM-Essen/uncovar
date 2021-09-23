rule simulate_strain_reads:
    input:
        get_genome_fasta,
    output:
        left=temp("resources/benchmarking/{accession}/reads.1.fastq.gz"),
        right=temp("resources/benchmarking/{accession}/reads.2.fastq.gz"),
    params:
        no_reads=lambda wildcards: no_reads(wildcards),
    log:
        "logs/mason/benchmarking/{accession}.log",
    conda:
        "../envs/mason.yaml"
    shell:  # median reads in data: 584903
        "mason_simulator -ir {input} -n {params.no_reads} -o {output.left} -or {output.right} 2> {log}"


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
        call_plots=expand(
            "results/benchmarking/plots/strain-calls/non-cov2-{accession}.strains.{caller}.svg",
            accession=get_non_cov2_accessions(),
            caller=["pangolin", "kallisto"],
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
        "results/benchmarking/plots/{caller}-strain-call-error-heatmap.svg",
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
        temp("results/{date}/assembly/{sample}/trinity/{sample}.contigs.fasta"),
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
        "mv {params.outdir}/Trinity.fasta {output} ) > {log} 2>&1"


rule assembly_comparison_velvet:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        temp("results/{date}/assembly/{sample}/velvet/{sample}.contigs.fasta"),
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
        contigs="results/{date}/assembly/{sample}/{assembler}/{sample}.contigs.fasta",
        reference="resources/genomes/main.fasta",
    output:
        temp(
            "results/{date}/assembly/{sample}/{assembler}/{sample}.ordered.contigs.fasta"
        ),
    log:
        "logs/{date}/ragoo/{assembler}/{sample}.log",
    params:
        outdir=get_output_dir,
    conda:
        "../envs/ragoo.yaml"
    shadow:
        "minimal"
    shell:
        "(cd {params.outdir} &&"
        " ragoo.py ../../../../../{input.contigs} ../../../../../{input.reference} &&"
        " cd ../../../../../ && mv {params.outdir}/ragoo_output/ragoo.fasta {output})"
        " > {log} 2>&1"


use rule filter_chr0 as filter_chr0_assembly_comparison with:
    input:
        "results/{date}/assembly/{sample}/{assembler}/{sample}.ordered.contigs.fasta",
    output:
        "results/{date}/assembly/{sample}/{assembler}/{sample}.contigs.ordered.filtered.fasta",
    log:
        "logs/{date}/ragoo/{assembler}/{sample}_cleaned.log",


use rule align_contigs as align_contigs_assembly_comparison with:
    input:
        target="resources/genomes/main.fasta",
        query="results/{date}/assembly/{sample}/{assembler}/{sample}.contigs.fasta",
    output:
        "results/{date}/assembly/{sample}/{assembler}/main_{sample}.bam",
    log:
        "results/{date}/assembly/{sample}/{assembler}/main_{sample}.log",


use rule quast as quast_assembly_comparison with:
    input:
        fasta="results/{date}/assembly/{sample}/{assembler}/{sample}.contigs.fasta",
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
            "results/{zip1}/assembly/{zip2}/{{exp}}/{zip2}.contigs.fasta"
        ),
        final=get_samples_for_assembler_comparison(
            "results/{zip1}/assembly/{zip2}/{{exp}}/{zip2}.contigs.ordered.filtered.fasta",
        ),
        quast=get_samples_for_assembler_comparison(
            "results/{zip1}/assembly/{zip2}/{{exp}}/quast/transposed_report.tsv",
        ),
    output:
        "results/benchmarking/plots/assembler-comparison.svg",
        "results/benchmarking/plots/assembler-comparison.csv",
        "results/benchmarking/plots/assembler-comparison_genome_fraction.svg",
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
        "results/benchmarking/plots/{caller}-call-dependency.svg",
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
