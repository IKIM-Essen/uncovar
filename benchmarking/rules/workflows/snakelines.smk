# source: https://snakelines.readthedocs.io/en/latest/pipelines/covidseq.html


rule snakelines_download:
    output:
        repo=directory("resources/benchmarking/snakelines/repo"),
        config="resources/benchmarking/snakelines/repo/example/covseq/covseq.yaml",
        snakeline="resources/benchmarking/snakelines/repo/snakelines.snake",
    log:
        "logs/snakelines_download.log",
    conda:
        "../envs/git.yaml"
    shell:
        "if [ -d '{output.repo}' ]; then rm -Rf {output.repo}; fi &&"
        "git clone https://github.com/thomasbtf/snakelines.git {output.repo} 2> {log}"


rule snakelines_prepare_primer:
    input:
        fasta="resources/genomes/{id}.fasta".format(
            id=config["preprocessing"]["amplicon-reference"]
        ),
        bed=config["preprocessing"]["amplicon-primers"],
    output:
        "resources/benchmarking/primer.fasta",
    log:
        "logs/snakelines_prepare_primer.log",
    conda:
        "../../envs/tools.yaml"
    shell:
        "bedtools getfasta  -fi {input.fasta} -bed {input.bed} -fo {output} 2> {log}"


rule snakelines_prepare_data:
    input:
        primers="resources/benchmarking/primer.fasta",
        samples=get_fastqs,
        reference="resources/genomes/main.fasta",
    output:
        outdir=temp(directory("results/benchmarking/snakelines/{sample}/")),
        fq1="results/benchmarking/snakelines/{sample}/reads/original/{sample}_R1.fastq.gz",
        fq2="results/benchmarking/snakelines/{sample}/reads/original/{sample}_R2.fastq.gz",
        reference="results/benchmarking/snakelines/{sample}/reference/sars_cov_2/sars_cov_2.fa",
        primers="results/benchmarking/snakelines/{sample}/reference/sars_cov_2/adapters/artic.fa",
    log:
        "logs/snakelines_prepare_data/{sample}.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "(ln -sr {input.samples[0]} {output.fq1} &&"
        " ln -sr {input.samples[1]} {output.fq2} &&"
        " ln -sr {input.reference} {output.reference} &&"
        " ln -sr {input.primers} {output.primers})"
        "2> {log}"


rule snakeline_unzip_hg38:
    input:
        "resources/genomes/human-genome.fna.gz",
    output:
        "resources/genomes/human-genome.fna",
    log:
        "logs/snakeline_unzip_hg38.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "gzip -dk {input} 2> {log}"


rule snakelines_prepare_hg38:
    input:
        "resources/genomes/human-genome.fna",
    output:
        "results/benchmarking/snakelines/{sample}/reference/hg38/hg38.fa",
    log:
        "logs/snakelines_prepare_hg38/{sample}.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "ln -sr {input} {output} 2> {log}"


rule snakeline:
    input:
        fq1="results/benchmarking/snakelines/{sample}/reads/original/{sample}_R1.fastq.gz",
        fq2="results/benchmarking/snakelines/{sample}/reads/original/{sample}_R1.fastq.gz",
        reference="results/benchmarking/snakelines/{sample}/reference/sars_cov_2/sars_cov_2.fa",
        primers="results/benchmarking/snakelines/{sample}/reference/sars_cov_2/adapters/artic.fa",
        config="resources/benchmarking/snakelines/repo/example/covseq/covseq.yaml",
        snakefile="resources/benchmarking/snakelines/repo/snakelines.snake",
    output:
        vcf="results/benchmarking/snakelines/{sample}/variant/sars_cov_2-wgs/original/{sample}.vcf",
        consensus="results/benchmarking/snakelines/{sample}/report/public/01-example/{sample}/consensus-sars_cov_2-wgs.fa",
        pangolin="results/benchmarking/snakelines/{sample}/report/public/01-example/{sample}/lineage_report-sars_cov_2-wgs.csv",
    log:
        "logs/snakeline/{sample}.log",
    conda:
        "../../envs/snakelines.yaml"
    benchmark:
        "benchmarks/snakeline/{sample}.benchmark.txt"
    threads: 4
    resources:
        snakeline=1,
    params:
        outdir=lambda w: f"results/benchmarking/snakelines/{w.sample}",
        cwd=lambda w: os.getcwd(),
    shell:
        "(cd {params.outdir} && "
        " snakemake --snakefile {params.cwd}/{input.snakefile}"
        " --cores {threads}"
        " --configfile {params.cwd}/{input.config}"
        " --verbose"
        " --config threads={threads}"
        " --rerun-incomplete"
        " --use-conda)"
        "> {log} 2>&1"
