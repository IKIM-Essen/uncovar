# source: https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe#3-usage
rule CovPipe_prepare_samples:
    input:
        get_fastqs,
    output:
        directory("resources/benchmarking/data/CovPipe/{sample}"),
    log:
        "logs/CovPipe_prepare_samples/{sample}.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "mkdir {output} && cp {input[0]} {output} && cp {input[1]} {output}"


rule CovPipe_prepare_adapter_file:
    output:
        "resources/benchmarking/data/CovPipe/adapters/{sample}/adapters.fasta",
    log:
        "logs/CovPipe_prepare_adapter_file/{sample}.log",
    conda:
        "../../envs/unix.yaml"
    params:
        adapters=lambda w: get_adapters(w)
        .replace("--adapter_sequence ", ">adapter fwd\n")
        .replace(" --adapter_sequence_r2 ", "\n>adapter rev\n"),
    shell:
        "echo '{params.adapters}' > {output}"


rule CovPipe:
    input:
        input_dir="resources/benchmarking/data/CovPipe/{sample}/",
        reference="resources/genomes/main.fasta",
        adapter="resources/benchmarking/data/CovPipe/adapters/{sample}/adapters.fasta",
        primer="resources/primer.bedpe",
    output:
        out_dir=directory("results/benchmarking/CovPipe/{sample}-{covpipe_name}/"),
        consensuses_masked="results/benchmarking/CovPipe/{sample}-{covpipe_name}/results/consensuses_masked/{covpipe_name}.masked_consensus.fasta",
        consensuses_iupac="results/benchmarking/CovPipe/{sample}-{covpipe_name}/results/consensuses_iupac/{covpipe_name}.iupac_consensus.fasta",
        vcf="results/benchmarking/CovPipe/{sample}-{covpipe_name}/results/intermediate_data/04_variant_calling/{covpipe_name}/{covpipe_name}.vcf",
        pangolin="results/benchmarking/CovPipe/{sample}-{covpipe_name}/results/intermediate_data/06_lineages/all_samples.lineage.txt",
    resources:
        external_pipeline=1,
    log:
        "logs/CovPipe/{sample}-{covpipe_name}.log",
    conda:
        "../../envs/covpipe.yaml"
    shell:
        "(ncov_minipipe"
        " --reference {input.reference}"
        " --input {input.input_dir}"
        " --adapter {input.adapter}"
        " --primer {input.primer}"
        " --pangolin $(which pangolin | sed 's/\/bin.*//g')"
        " -o {output.out_dir})"
        " > {log} 2>&1"


rule CovPipe_results:
    input:
        lambda w: expand(
            [
                "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/consensuses_masked/{covpipe_name}.masked_consensus.fasta",
                "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/consensuses_iupac/{covpipe_name}.iupac_consensus.fasta",
                "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/intermediate_data/04_variant_calling/{covpipe_name}/{covpipe_name}.vcf",
                "results/benchmarking/CovPipe/{{sample}}-{covpipe_name}/results/intermediate_data/06_lineages/all_samples.lineage.txt",
            ],
            covpipe_name=w.sample.replace("_", "__"),
        ),
    output:
        touch("results/benchmarking/CovPipe/{sample}/.done"),
    log:
        "logs/CovPipe_results/{sample}.log",
