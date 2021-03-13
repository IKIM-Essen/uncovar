rule simulate_strain_reads:
    input:
        get_genome_fasta
    output:
        left="resources/benchmarking/{accession}/reads.1.fastq.gz",
        right="resources/benchmarking/{accession}/reads.2.fastq.gz",
    params:
        no_reads = lambda wildcards: no_reads(wildcards)
    log:
        "logs/mason/benchmarking/{accession}.log",
    conda:
        "../envs/mason.yaml"
    shell:  # median reads in data: 584903
        "mason_simulator -ir {input} -n {params.no_reads} -o {output.left} -or {output.right} 2> {log}"


rule mix_strain_reads:
    input:
        left = expand("resources/benchmarking/{mix}/reads.1.fastq.gz", mix = ["#{{strain_{}}}".format(i) for i in range(config["mixtures"]["no_strains"])]),
        right = expand("resources/benchmarking/{mix}/reads.2.fastq.gz", mix = ["#{{strain_{}}}".format(i) for i in range(config["mixtures"]["no_strains"])]),
    output:
        left = expand("resources/mixtures/{mix}/reads.1.fastq.gz", mix = "".join(["#{{strain_{}}}".format(i) for i in range(config["mixtures"]["no_strains"])])),
        right = expand("resources/mixtures/{mix}/reads.2.fastq.gz", mix = "".join(["#{{strain_{}}}".format(i) for i in range(config["mixtures"]["no_strains"])]))
    log:
        "logs/mix_strain_reads/{}".format("".join(["#{{strain_{}}}".format(i) for i in range(config["mixtures"]["no_strains"])]))
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
        "results/benchmarking/polished-contigs/benchmark-sample-{accession}.fasta",
    output:
        "results/benchmarking/assembly/{accession}.bam",
    log:
        "logs/test-assembly-results/{accession}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 --MD --eqx -ax asm5 {input} -o {output} 2> {log}"


rule summarize_assembly_results:
    input:
        bams=get_assembly_comparisons(bams=True),
        refs=get_assembly_comparisons(bams=False),
    output:
        "results/benchmarking/assembly.csv",
    log:
        "logs/assembly/assembly-results.log",
    conda:
        "../envs/pysam.yaml"
    notebook:
        "../notebooks/assembly-benchmark-results.py.ipynb"


rule test_non_cov2:
    input:
        pangolin=get_non_cov2_calls(from_caller="pangolin"),
        kallisto=get_non_cov2_calls(from_caller="kallisto"),
    output:
        "results/test-cases/non-sars-cov-2.csv",
    params:
        accessions=get_non_cov2_accessions(),
    log:
        "../logs/test-cases/summarize_non_cov2.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize-non-cov2.py"


rule report_non_cov2:
    input:
        summary="results/test-cases/non-sars-cov-2.csv",
        call_plots=expand(
            "results/test-cases/plots/strain-calls/non-cov2-{accession}.strains.{caller}.svg",
            accession=get_non_cov2_accessions(),
            caller=["pangolin", "kallisto"],
        ),
    output:
        report(
            directory("results/test-cases/html"),
            htmlindex="index.html",
            category="Test results",
        ),
    log:
        "../logs/report_non_cov2.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report -s '\t' {input.summary} {output}"


rule evaluate_kallistos_read_error:
    input:
       get_mixture_results
    output:
        "results/benchmarking/kallisto-read-error.csv",
    params:
     max_reads = config["mixtures"]["max_reads"]
    log:
        "logs/evaluate-kallistos-read-error.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/evaluate-kallistos-read-error.py"


rule plot_kallistos_read_error:
    input:
        "results/benchmarking/kallisto-read-error.csv",
    output:
        "results/benchmarking/kallisto-read-error.svg",
    log:
        "logs/plot-kallistos-read-error.log"
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-kallistos-read-error.py"