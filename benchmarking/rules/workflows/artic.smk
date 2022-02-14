# source: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
rule download_artic_primer_schemes:
    output:
        directory("resources/benchmarking/artic/repo"),
    log:
        "logs/download_artic_primer_schemes.log",
    conda:
        "../envs/git.yaml"
    shell:
        "git clone https://github.com/artic-network/artic-ncov2019.git {output} 2> {log}"


rule artic_guppyplex:
    input:
        get_fastq_pass_path_barcode,
    output:
        outdir=temp(directory("results/benchmarking/artic/guppyplex/{sample}/")),
        fasta="results/benchmarking/artic/guppyplex/{sample}/{sample}.fasta",
    log:
        "logs/artic_guppyplex/{sample}.log",
    conda:
        "../../envs/artic.yaml"
    benchmark:
        "benchmarks/artic_guppyplex/{sample}.benchmark.txt"
    resources:
        external_pipeline=1,
    shell:
        "artic guppyplex --min-length 400 --max-length 700"
        " --directory {input} --output {output.fasta}"
        " > {log} 2>&1"


rule artic_minion_nanopolish:
    input:
        sequencing_summary=lambda wildcards: get_seq_summary(wildcards),
        fast5=lambda wildcards: get_fast5_pass_path_barcode(wildcards),
        fasta="results/benchmarking/artic/guppyplex/{sample}/{sample}.fasta",
        repo="resources/benchmarking/artic/repo",
        guppy_dir="results/benchmarking/artic/guppyplex/{sample}/",
    output:
        outdir=temp(directory("results/benchmarking/artic/minion/nanopolish/{sample}/")),
        vcf="results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.merged.vcf",
        consensus="results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.consensus.fasta",
    log:
        "logs/artic_minion/nanopolish/{sample}.log",
    conda:
        "../../envs/artic.yaml"
    benchmark:
        "benchmarks/artic_nanopolish/{sample}.benchmark.txt"
    threads: 4
    params:
        primer_schemes=lambda w, input: os.path.join(input.repo, "primer_schemes"),
        cwd=lambda w: os.getcwd(),
    resources:
        external_pipeline=1,
    shell:
        "(cd {output.outdir} &&"
        " artic minion --normalise 200 --threads {threads}"
        " --scheme-directory {params.cwd}/{params.primer_schemes}"
        " --read-file {params.cwd}/{input.fasta}"
        " --fast5-directory {params.cwd}/{input.fast5}"
        " --sequencing-summary {params.cwd}/{input.sequencing_summary}"
        " nCoV-2019/V3 {wildcards.sample})"
        " > {params.cwd}/{log} 2>&1"


rule artic_minion_medaka:
    input:
        fast5=lambda wildcards: get_fast5_pass_path_barcode(wildcards),
        fasta="results/benchmarking/artic/guppyplex/{sample}/{sample}.fasta",
        repo="resources/benchmarking/artic/repo",
        guppy_dir="results/benchmarking/artic/guppyplex/{sample}/",
    output:
        outdir=temp(directory("results/benchmarking/artic/minion/medaka/{sample}/")),
        vcf="results/benchmarking/artic/minion/medaka/{sample}/{sample}.merged.vcf",
        consensus="results/benchmarking/artic/minion/medaka/{sample}/{sample}.consensus.fasta",
    log:
        "logs/artic_minion/medaka/{sample}.log",
    conda:
        "../../envs/artic.yaml"
    benchmark:
        "benchmarks/artic_medaka/{sample}.benchmark.txt"
    threads: 4
    params:
        primer_schemes=lambda w, input: os.path.join(input.repo, "primer_schemes"),
        medaka_model=config["assembly"]["oxford nanopore"]["medaka_model"],
        cwd=lambda w: os.getcwd(),
    resources:
        external_pipeline=1,
    shell:
        "(cd {output.outdir} &&"
        " artic minion --medaka --medaka-model {params.medaka_model}"
        " --normalise 200 --threads {threads}"
        " --scheme-directory {params.cwd}/{params.primer_schemes}"
        " --read-file {params.cwd}/{input.fasta}"
        " nCoV-2019/V3 {wildcards.sample})"
        " > {params.cwd}/{log} 2>&1"
