# rule uncovar_deploy_workflow:
#     output:
#         directory("results/benchmarking/UnCoVar/"),
#     log:
#         "logs/uncovar_deploy_workflow/.log"
#     conda:
#         "../../envs/snakedeploy.yaml"
#     shell:
#         "snakedeploy deploy-workflow https://github.com/IKIM-Essen/uncovar {output} --tag v0.14.0 > {log} 2&>1"


rule uncovar_download:
    output:
        touch("results/benchmarking/.indicators/uncovar/dowloaded"),
    log:
        "logs/uncovar_download.log",
    conda:
        "../../envs/git.yaml"
    params:
        repo=directory("results/benchmarking/UnCoVar/"),
    shell:
        "if [ -d '{params.repo}' ]; then rm -Rf {params.repo}; fi &&"
        "git clone --branch v0.14.0 https://github.com/IKIM-Essen/uncovar {params.repo} 2> {log}"


rule uncovar_prep_sample_sheet:
    input:
        "results/benchmarking/.indicators/uncovar/dowloaded",
    output:
        touch("results/benchmarking/.indicators/uncovar/samples_preped"),
    log:
        "logs/uncovar_prep_sample_sheet.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "cp benchmarking/config/samples/samples.csv results/benchmarking/UnCoVar/config/pep/samples.csv 2> {log}"


rule link_data:
    input:
        data="data/",
        uncvoar="results/benchmarking/.indicators/uncovar/dowloaded",
    output:
        touch("results/benchmarking/.indicators/uncovar/data_preped"),
    log:
        "logs/link_data.log",
    conda:
        "../../envs/unix.yaml"
    params:
        data_dir="results/benchmarking/UnCoVar/data",
    shell:
        "ln -sr {input.data} {params.data_dir}"


rule uncovar_run:
    input:
        uncovar_dir="results/benchmarking/.indicators/uncovar/dowloaded",
        samples="results/benchmarking/.indicators/uncovar/samples_preped",
        data="results/benchmarking/.indicators/uncovar/data_preped",
    output:
        "results/benchmarking/UnCoVar/results/{date}/filtered-calls/ref~main/{sample}.subclonal.high+moderate-impact.orf.bcf",
    log:
        "logs/uncovar_run/{sample}-{date}.log",
    conda:
        "../../envs/snakemake.yaml"
    benchmark:
        "benchmarks/uncovar/{sample}~{date}.benchmark.txt"
    threads: 4
    resources:
        uncovar=1,
    params:
        uncovar_dir="results/benchmarking/UnCoVar/",
        output=lambda w, output: " ".join(
            path.removeprefix("results/benchmarking/UnCoVar/") for path in output
        ),
    shell:
        "(cd {params.uncovar_dir} && "
        " snakemake --cores {threads} --use-conda --resources ncbi_api_requests=1 --rerun-incomplete {params.output}"
        ") > {log} 2>&1"


rule uncovar_bcf_2_vcf:
    input:
        lambda w: expand(
            "results/benchmarking/UnCoVar/results/{date}/filtered-calls/ref~main/{{sample}}.subclonal.high+moderate-impact.orf.bcf",
            date=get_date_for_sample(w),
        ),
    output:
        "results/benchmarking/UnCoVar/benchmark-result/{sample}.vcf",
    log:
        "logs/uncovar_bcf_2_vcf/{sample}.log",
    conda:
        "../../envs/tools.yaml"
    shell:
        "bcftools view {input} > {output}"


# TODO: check if the other pipelines filter for impact or other criteria
