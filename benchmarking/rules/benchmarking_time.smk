rule annotate_time_benchmark:
    input:
        "results/benchmarking/backups/{output_type}/{tech}/{workflow}/{sample}.some.extension",
    output:
        "results/benchmarking/tabels/time_benchmarks/{output_type}/{tech}/{workflow}/{sample}.tsv",
    log:
        "logs/annotate_time_benchmark/{output_type}/{tech}/{workflow}/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/annotate_time_output.py"


rule agg_time_benchmark:
    input:
        get_output_of_pipelines(
            path="results/benchmarking/tabels/time_benchmarks/{output_type}/{tech}/{workflow}/{sample}.tsv",
            output="time",
        ),
    output:
        "results/benchmarking/tabels/time_benchmarks/time.csv",
    log:
        "logs/agg_time_benchmark.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/agg_workflow_output.py"


rule get_size_of_input:
    input:
        get_paths_of_input_files,
    output:
        "results/benchmarking/tabels/time_benchmarks/input_size/{sample}.tsv",
    log:
        "logs/get_size_of_input/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/count_basepairs.py"


rule plot_time:
    input:
        times="results/benchmarking/tabels/time_benchmarks/time.csv",
        sizes=expand(
            "results/benchmarking/tabels/time_benchmarks/input_size/{sample}.tsv",
            sample=get_samples(),
        ),
    output:
        "results/benchmarking/plots/time.svg",
    log:
        "logs/plot_time.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_time.py"
