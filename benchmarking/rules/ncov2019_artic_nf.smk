# source: https://github.com/connor-lab/ncov2019-artic-nf
rule ncov2019_artic_nf_illumina_data_prep:
    input:
        get_fastqs,
    output:
        d=temp(directory("resources/data/ncov2019-artic-nf/illumina/{sample}")),
        fq1=temp(
            "resources/data/ncov2019-artic-nf/illumina/{sample}/{sample}_R1.fastq.gz"
        ),
        fq2=temp(
            "resources/data/ncov2019-artic-nf/illumina/{sample}/{sample}_R2.fastq.gz"
        ),
    log:
        "logs/ncov2019_artic_nf_illumina_data_prep/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {output.d} &&"
        " cp {input[0]} {output.fq1} &&"
        " cp {input[1]} {output.fq2})"
        " 2> {log}"


rule ncov2019_artic_nf_illumina:
    input:
        directory="resources/data/ncov2019-artic-nf/illumina/{sample}/",
    output:
        outdir=temp(
            directory("results/benchmarking/ncov2019_artic_nf/illumina/{sample}/")
        ),
        consensus="results/benchmarking/ncov2019_artic_nf/illumina/{sample}/ncovIllumina_sequenceAnalysis_makeConsensus/{sample}.primertrimmed.consensus.fa",
        vcf="results/benchmarking/ncov2019_artic_nf/illumina/{sample}/ncovIllumina_sequenceAnalysis_callVariants/{sample}.variants.tsv",
    log:
        "logs/ncov2019_artic_nf/illumina/{sample}.log",
    threads: 8
    params:
        pipeline="connor-lab/ncov2019-artic-nf",
        revision="v1.3.0",
        profile=["conda"],
        flags="--illumina",
        outdir=lambda w: f"results/benchmarking/ncov2019_artic_nf/illumina/{w.sample}",
        prefix=lambda w: w.sample,
    handover: True
    conda:
        "../envs/nextflow.yaml"
    script:
        "../scripts/nextflow.py"


rule ncov2019_artic_nf_nanopore_data_prep:
    input:
        get_fastq_or_fast5,
    output:
        temp(
            directory(
                "resources/benchmarking/data/ncov2019_artic_nf/nanopore/{sample}/{folder}/"
            )
        ),
    log:
        "logs/ncov2019_artic_nf_nanopore_data_prep/{sample}-{folder}.log",
    conda:
        "../envs/unix.yaml"
    params:
        barcode=lambda w, output: os.path.join(output[0], get_barcode(w)),
    shell:
        "mkdir -p {params.barcode} && cp -r {input} {output}"


use rule ncov2019_artic_nf_illumina as ncov2019_artic_nf_nanopore_nanopolish with:
    input:
        basecalled_fastq="resources/benchmarking/data/ncov2019_artic_nf/nanopore/{sample}/fastq_pass/",
        fast5_pass="resources/benchmarking/data/ncov2019_artic_nf/nanopore/{sample}/fast5_pass/",
        sequencing_summary=lambda wildcards: get_seq_summary(wildcards),
    output:
        outdir=temp(
            directory(
                "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{sample}/"
            )
        ),
        consensus=temp(
            "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{sample}/articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish/{sample}_{barcode}.consensus.fasta"
        ),
        vcf=temp(
            "results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{sample}/articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish/{sample}_{barcode}.merged.vcf"
        ),
    log:
        "logs/ncov2019_artic_nf/nanopore/nanopolish/{sample}-{barcode}.log",
    params:
        pipeline="connor-lab/ncov2019-artic-nf",
        revision="v1.3.0",
        flags="--nanopolish",
        outdir=lambda w: f"results/benchmarking/ncov2019_artic_nf/nanopore/nanopolish/{w.sample}",
        prefix=lambda w: w.sample,
    conda:
        "../envs/nextflow_ncov2019_artic_nf_nanopore.yaml"


use rule ncov2019_artic_nf_nanopore_nanopolish as bncov2019_artic_nf_nanopore_medaka with:
    output:
        outdir=temp(
            directory(
                "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{sample}/"
            )
        ),
        consensus=temp(
            "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{sample}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/{sample}_{barcode}.consensus.fasta"
        ),
        vcf=temp(
            "results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{sample}/articNcovNanopore_sequenceAnalysisMedaka_articMinIONMedaka/{sample}_{barcode}.consensus.fasta"
        ),
    log:
        "logs/ncov2019_artic_nf/nanopore/medaka/{sample}-{barcode}.log",
    params:
        pipeline="connor-lab/ncov2019-artic-nf",
        revision="v1.3.0",
        flags="--medaka",
        outdir=lambda w: f"results/benchmarking/ncov2019_artic_nf/nanopore/medaka/{w.sample}",
        prefix=lambda w: w.sample,
