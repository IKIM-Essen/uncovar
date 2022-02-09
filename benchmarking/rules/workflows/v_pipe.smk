# source: https://cbg-ethz.github.io/V-pipe/tutorial/sars-cov2/
rule v_pipe_work_dir:
    input:
        "resources/benchmarking/v-pipe/repo",
    output:
        touch("results/benchmarking/v-pipe/{sample}/.work-dir-created"),
    log:
        "logs/v_pipe_work_dir/{sample}.log",
    conda:
        "../../envs/v-pipe.yaml"
    params:
        v_pipe_copy="results/benchmarking/v-pipe/{sample}/",
        work_dir_path="results/benchmarking/v-pipe/{sample}/work",
    shell:
        "(ln -sr {input}/* {params.v_pipe_copy} &&"
        " mkdir -p {params.work_dir_path} &&"
        " cd {params.work_dir_path} &&"
        " ../init_project.sh)"
        " > {log} 2>&1"


rule v_pipe_setup_samples:
    input:
        workdir="results/benchmarking/v-pipe/{sample}/.work-dir-created",
        fastqs=get_fastqs,
    output:
        expand(
            "results/benchmarking/v-pipe/{{sample}}/work/samples/{{sample}}/20200102/raw_data/{{sample}}_R{read}.fastq",
            read=[1, 2],
        ),
    log:
        "logs/v_pipe_setup_samples/{sample}.log",
    conda:
        "../../envs/v-pipe.yaml"
    params:
        fq_dir=lambda w, input: os.path.join(
            os.path.dirname(input.workdir),
            "work",
            "samples",
            w.sample,
            "20200102",
            "raw_data",
        ),
    shell:
        "(mkdir -p {params.fq_dir} &&"
        " gzip -d {input.fastqs[0]} -c > {output[0]} &&"
        " gzip -d {input.fastqs[1]} -c > {output[1]})"
        " 2> {log}"


rule v_pipe_dry_run:
    input:
        expand(
            "results/benchmarking/v-pipe/{{sample}}/work/samples/{{sample}}/20200102/raw_data/{{sample}}_R{read}.fastq",
            read=[1, 2],
        ),
    output:
        "results/benchmarking/v-pipe/{sample}/work/samples.tsv",
    log:
        "logs/v_pipe_dry_run/{sample}.log",
    conda:
        "../../envs/v-pipe.yaml"
    params:
        workdir=lambda w: f"results/benchmarking/v-pipe/{w.sample}/work",
    resources:
        external_pipeline=1,
    shell:
        "(cd {params.workdir} &&"
        " ./vpipe --dryrun)"
        " > {log} 2>&1"


rule v_pipe_update_sample_sheet:
    input:
        "results/benchmarking/v-pipe/{sample}/work/samples.tsv",
    output:
        touch("results/benchmarking/v-pipe/{sample}/.edited-sample"),
    log:
        "logs/v_pipe_update_sample_sheet/{sample}.log",
    conda:
        "../../envs/v-pipe.yaml"
    shell:
        "sed -i 's/$/\t150/' {input} 2> {log}"


rule v_pipe_run:
    input:
        "results/benchmarking/v-pipe/{sample}/.edited-sample",
    output:
        vcf="results/benchmarking/v-pipe/{sample}/work/samples/{sample}/20200102/variants/SNVs/snvs.vcf",
        consensus="results/benchmarking/v-pipe/{sample}/work/samples/{sample}/20200102/references/ref_majority.fasta",
    log:
        "logs/v_pipe_run/{sample}.log",
    conda:
        "../../envs/v-pipe.yaml"
    benchmark:
        "benchmarks/v_pipe/{sample}.benchmark.txt"
    threads: 16
    resources:
        external_pipeline=1,
    params:
        workdir=lambda w, input: os.path.join(os.path.dirname(input[0]), "work"),
    shell:
        "(cd {params.workdir} &&"
        " ./vpipe --cores {threads} -p -F)"
        " > {log} 2>&1"


rule v_pipe_fix_vcf:
    input:
        "results/benchmarking/v-pipe/{sample}/work/samples/{sample}/20200102/variants/SNVs/snvs.vcf",
    output:
        "results/benchmarking/v-pipe/{sample}/work/samples/{sample}/20200102/variants/SNVs/fixed-vcf/snvs.vcf",
    log:
        "logs/v_pipe_fix_vcf/{sample}.log",
    conda:
        "../../envs/python.yaml"
    script:
        "../../scripts/v_pipe_fix_vcf.py"
