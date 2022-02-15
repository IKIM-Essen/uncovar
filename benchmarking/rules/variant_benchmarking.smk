rule rename_contig:
    input:
        vcf=get_output_from_pipline("vcf"),
        outdir=get_output_from_pipline("outdir"),
    output:
        "results/benchmarking/variant-calls/renamed/{workflow}/{sample}.vcf",
    log:
        "logs/rename_contig/{workflow}/{sample}.log",
    params:
        search_string="MN908947.3",
        replace_string="NC_045512.2",
    conda:
        "../envs/unix.yaml"
    shell:
        "sed 's/{params.search_string}/{params.replace_string}/' {input.vcf} > {output} 2> {log}"


rule check_contig_flag:
    input:
        "results/benchmarking/variant-calls/renamed/{workflow}/{sample}.vcf",
    output:
        "results/benchmarking/variant-calls/contig-checked/{workflow}/{sample}.vcf",
    log:
        "logs/check_contig_flag/{workflow}/{sample}.log",
    params:
        search_string="##contig=<ID=NC_045512.2>",
        add_string="##contig=<ID=NC_045512.2>",
    conda:
        "../envs/unix.yaml"
    shell:
        "if ! grep -Fxq '{params.search_string}' {input}; "
        "then sed '2i {params.add_string}' {input} > {output}; "
        "else cp {input} {output}; fi"


rule check_genotype:
    input:
        "results/benchmarking/variant-calls/contig-checked/{workflow}/{sample}.vcf",
    output:
        "results/benchmarking/variant-calls/genotyped/{workflow}/{sample}.vcf",
    log:
        "logs/check_genotype/{workflow}/{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/check_genotype.py"


rule normalize_calls:
    input:
        "results/benchmarking/variant-calls/genotyped/{workflow}/{sample}.vcf",
        genome="resources/genomes/main.fasta",
        genome_index="resources/genomes/main.fasta.fai",
    output:
        "results/benchmarking/variant-calls/normalized-variants/{workflow}/{sample}.vcf.gz",
    params:
        extra=lambda w, input: f"--atomize -f {input.genome} --rm-dup exact",
    log:
        "logs/normalize-calls/{workflow}/{sample}.log",
    conda:
        "../envs/tools.yaml"
    wrapper:
        "v1.0.0/bio/bcftools/norm"


rule stratify:
    input:
        variants="results/benchmarking/variant-calls/normalized-variants/{workflow}/{sample}.vcf.gz",
        regions="results/benchmarking/sanger/aligned/{sample}.bed",  # we only have sanger as truth
    output:
        "results/benchmarking/variant-calls/stratified/{workflow}/{sample}.vcf.gz",
    log:
        "logs/stratify-truth/{workflow}/{sample}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "bedtools intersect -b {input.regions} -a <(bcftools view {input.variants}) -wa -f 1.0 -header | bcftools view -Oz > {output} 2> {log}"


rule bcftools_index:
    input:
        "results/benchmarking/{infix}.vcf.gz",
    output:
        "results/benchmarking/{infix}.vcf.gz.csi",
    log:
        "logs/bcftools-index/{infix}.log",
    wrapper:
        "v1.0.0/bio/bcftools/index"


rule benchmark_variants:
    input:
        truth="results/benchmarking/variant-calls/{source}/{workflow_A}/{sample}.vcf.gz",
        truth_idx="results/benchmarking/variant-calls/{source}/{workflow_A}/{sample}.vcf.gz.csi",
        query="results/benchmarking/variant-calls/{source}/{workflow_B}/{sample}.vcf.gz",
        genome="resources/genomes/main.fasta",
        genome_index="resources/genomes/main.fasta.fai",
    output:
        multiext(
            "results/benchmarking/happy/{source}/{workflow_A}-vs-{workflow_B}/{sample}/report",
            ".runinfo.json",
            ".vcf.gz",
            ".summary.csv",
            ".extended.csv",
            ".metrics.json.gz",
            ".roc.all.csv.gz",
            ".roc.Locations.SNP.csv.gz",
        ),
        # ".roc.Locations.INDEL.csv.gz", # bc of empty vcfs files
        # ".roc.Locations.INDEL.PASS.csv.gz",
    params:
        prefix=lambda w, input, output: output[0].split(".")[0],
        engine="vcfeval",
    log:
        "logs/happy/{source}/{workflow_A}-vs-{workflow_B}/{sample}.log",
    wrapper:
        "v1.0.0/bio/hap.py/hap.py"


rule agg_happy:
    input:
        get_benchmark_path(
            "results/benchmarking/happy/{{source}}/sanger-vs-{workflow}/{sample}/report.summary.csv"
        ),
    output:
        "results/benchmarking/workflow-comparison.{source}.tsv",
    log:
        "logs/agg_happy/{source}.log",
    conda:
        "../envs/python.yaml"
    params:
        metadata=get_benchmark_path("{workflow},{sample}"),
        platforms=get_benchmark_platforms,
    script:
        "../scripts/agg_happy.py"


rule agg:
    input:
        expand(
            "results/benchmarking/workflow-comparison.{source}.tsv",
            source=["stratified", "normalized-variants"],
        ),


rule plot_plot_precision_recall:
    input:
        "results/benchmarking/workflow-comparison.stratified.tsv",
    output:
        plot="results/benchmarking/workflow-comparison.svg",
    log:
        "logs/plot_plot_precision_recall.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_plot_precision_recall.py"
