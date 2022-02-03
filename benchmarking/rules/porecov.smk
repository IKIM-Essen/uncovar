# source: https://github.com/replikation/poreCov
rule poreCov_sample_sheet:
    input:
        get_fastq_pass_path_barcode,
    output:
        "results/benchmarking/poreCov/sample-sheets/{sample}/sample_names.csv",
    log:
        "logs/poreCov_sample_sheet/{sample}.log",
    conda:
        "../envs/unix.yaml"
    params:
        barcode=lambda w, input: os.path.basename(input[0]),
    shell:
        "(echo _id,Status,Description > {output} &&"
        " echo {wildcards.sample},{params.barcode},- >> {output})"
        " 2> {log}"


rule poreCov:
    input:
        fastq_pass=get_fastq_pass_path_barcode,
        sample_names="results/benchmarking/poreCov/sample-sheets/{sample}/sample_names.csv",
    output:
        outdir=directory("results/benchmarking/poreCov/{sample}/"),
        consensus="results/benchmarking/poreCov/{sample}/2.Genomes/all_consensus_sequences/{sample}.consensus.fasta",
        lineage_call="results/benchmarking/poreCov/{sample}/3.Lineages_Clades_Mutations/{sample}/lineage_report_{sample}.csv",
        # varaints=""
    log:
        "logs/poreCov/{sample}.log",
    params:
        pipeline="replikation/poreCov",
        revision="1.0.0",
        qs=lambda w, threads: threads,
        profile=["local", "docker"],
        flags="--update",
        cores=lambda w, threads: threads,
        output=lambda w: f"results/benchmarking/poreCov/{w.sample}",
        samples=lambda W, input: input.sample_names,
    handover: True
    threads: 16
    conda:
        "../envs/nextflow.yaml"
    resources:
        external_pipeline=1,
        nextflow=1,
    script:
        "../scripts/nextflow.py"
