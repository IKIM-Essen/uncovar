rule annotate_pangolin_calls:
    input:
        "results/benchmarking/backups/{output_type}/{tech}/{workflow}/{sample}.some.extension",
    output:
        "results/benchmarking/tabels/{output_type}/{tech}/{workflow}/{sample}.csv",
    log:
        "logs/annotate_pangolin_calls/{output_type}/{tech}/{workflow}/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/annotate_pangolin_output.py"


rule agg_pangolin_calls:
    input:
        get_output_of_pipelines(
            path="results/benchmarking/tabels/{output_type}/{tech}/{workflow}/{sample}.csv",
            output="pangolin",
        ),
    output:
        "results/benchmarking/tabels/pangolin.csv",
    log:
        "logs/agg_pangolin_calls.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/agg_workflow_output.py"


rule plot_pangolin_calls:
    input:
        "results/benchmarking/tabels/pangolin.csv",
    output:
        "results/benchmarking/plots/pangolin.svg",
    log:
        "logs/plot_pangolin_calls.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_pangolin_calls.py"
