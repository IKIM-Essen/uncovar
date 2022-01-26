# source: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
rule artic_guppyplex:
    input:
        get_fastq_pass_path_barcode,
    output:
        "results/benchmarking/artic/guppyplex/{sample}.fasta",
    log:
        "logs/artic_guppyplex/{sample}.log",
    conda:
        "../envs/artic.yaml"
    shell:
        "artic guppyplex --min-length 400 --max-length 700"
        " --directory {input} --output {output}"
        " > {log} 2>&1"


rule artic_minion_nanopolish:
    input:
        sequencing_summary=lambda wildcards: get_seq_summary(wildcards),
        fast5=lambda wildcards: get_fast5_pass_path_barcode(wildcards),
        fasta="results/benchmarking/artic/guppyplex/{sample}.fasta",
        repo="resources/benchmarking/artic/repo",
    output:
        vcf="results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.merged.vcf",
        consensus="results/benchmarking/artic/minion/nanopolish/{sample}/{sample}.consensus.fasta",
    log:
        "logs/artic_minion/nanopolish/{sample}.log",
    threads: 16
    conda:
        "../envs/artic.yaml"
    params:
        primer_schemes=lambda w, input: os.path.join(input.repo, "primer_schemes"),
        medaka_model=config["assembly"]["oxford nanopore"]["medaka_model"],
        outdir=get_output_dir,
        cwd=lambda w: os.getcwd(),
    shell:
        "(cd {params.outdir} &&"
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
        fasta="results/benchmarking/artic/guppyplex/{sample}.fasta",
        repo="resources/benchmarking/artic/repo",
    output:
        vcf="results/benchmarking/artic/minion/medaka/{sample}/{sample}.merged.vcf",
        consensus="results/benchmarking/artic/minion/medaka/{sample}/{sample}.consensus.fasta",
    log:
        "logs/artic_minion/medaka/{sample}.log",
    threads: 16
    conda:
        "../envs/artic.yaml"
    params:
        primer_schemes=lambda w, input: os.path.join(input.repo, "primer_schemes"),
        medaka_model=config["assembly"]["oxford nanopore"]["medaka_model"],
        outdir=get_output_dir,
        cwd=lambda w: os.getcwd(),
    shell:
        "(cd {params.outdir} &&"
        " artic minion --medaka --medaka-model {params.medaka_model}"
        " --normalise 200 --threads {threads}"
        " --scheme-directory {params.cwd}/{params.primer_schemes}"
        " --read-file {params.cwd}/{input.fasta}"
        " nCoV-2019/V3 {wildcards.sample})"
        " > {params.cwd}/{log} 2>&1"
