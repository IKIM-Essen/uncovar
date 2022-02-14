# source: https://bitbucket.org/auto_cov_pipeline/havoc/src/master/


rule HaVoc_prepare_ref:
    input:
        reference="resources/genomes/main.fasta",
    output:
        "resources/benchmarking/havoc/ref.fa",
    log:
        "logs/HaVoc_prepare_ref.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "ln -sr {input.reference} {output[0]} 2>{log}"


rule HaVoc_prepare_adapter_file:
    output:
        "resources/benchmarking/havoc/NexteraPE-PE.fa",
    log:
        "logs/HaVoc_prepare_adapter_file.log",
    conda:
        "../../envs/unix.yaml"
    params:
        adapters=lambda w: get_adapters(w)
        .replace("--adapter_sequence ", ">adapter fwd\n")
        .replace(" --adapter_sequence_r2 ", "\n>adapter rev\n"),
    shell:
        "echo '{params.adapters}' > {output} 2> {log}"


rule HaVoc_prepare_data:
    input:
        get_fastqs,
    output:
        fq1="results/benchmarking/havoc/{sample}/data/{sample}_R1.fastq.gz",
        fq2="results/benchmarking/havoc/{sample}/data/{sample}_R2.fastq.gz",
    log:
        "logs/HaVoc_prepare_data/{sample}.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "(ln -sr {input[0]} {output.fq1} &&"
        " ln -sr {input[1]} {output.fq2})"
        "2>{log}"


rule HaVoc:
    input:
        script="resources/benchmarking/havoc/HAVoC.sh",
        fq1="results/benchmarking/havoc/{sample}/data/{sample}_R1.fastq.gz",
        fq2="results/benchmarking/havoc/{sample}/data/{sample}_R2.fastq.gz",
        ref="resources/benchmarking/havoc/ref.fa",
        adapters="resources/benchmarking/havoc/NexteraPE-PE.fa",
    output:
        out_dir=temp(
            directory("results/benchmarking/havoc/{sample}/data/{havoc_name}/")
        ),
        consensus="results/benchmarking/havoc/{sample}/data/{havoc_name}/{havoc_name}_consensus.fa",
        pangolin="results/benchmarking/havoc/{sample}/data/{havoc_name}/{havoc_name}_pangolin_lineage.csv",
        vcf="results/benchmarking/havoc/{sample}/data/{havoc_name}/{havoc_name}_indel.vcf",
    log:
        "logs/HaVoc/{sample}-{havoc_name}.log",
    conda:
        "../../envs/havoc.yaml"
    benchmark:
        "benchmarks/havoc/{sample}~{havoc_name}.benchmark.txt"
    threads: 4
    params:
        fastq_dir=lambda w, input: os.path.dirname(input.fq1),
        out_dir=lambda w, input: os.path.dirname(os.path.dirname(input.fq1)),
        cwd=os.getcwd(),
    shell:
        "(cd {params.out_dir} &&"
        " bash {params.cwd}/{input.script} {params.cwd}/{params.fastq_dir}) > {log} 2>&1"


rule HaVoc_fix_vcf:
    input:
        lambda w: "results/benchmarking/havoc/{{sample}}/data/{havoc_name}/{havoc_name}_indel.vcf".format(
            havoc_name=w.sample.split("_")[0]
        ),
    output:
        "results/benchmarking/havoc/{sample}/{sample}.fixed.vcf",
    log:
        "logs/HaVoc_fix_vcf/{sample}.log",
    params:
        search_string=lambda w: w.sample.split("_")[0],
        replace_string="NC_045512.2",
    conda:
        "../../envs/unix.yaml"
    shell:
        "sed 's/{params.search_string}/{params.replace_string}/' {input} > {output} 2> {log}"
