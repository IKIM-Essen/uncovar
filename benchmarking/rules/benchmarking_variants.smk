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
        vcf="results/benchmarking/variant-calls/genotyped/{workflow}/{sample}.vcf",
        genome="resources/genomes/main.fasta",
        genome_index="resources/genomes/main.fasta.fai",
    output:
        "results/benchmarking/variant-calls/normalized-variants/{workflow}/{sample}.vcf.gz",
    log:
        "logs/normalize-calls/{workflow}/{sample}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(bcftools norm --check-ref s  -f {input.genome} {input.vcf} |"
        " bcftools norm --output-type z -o {output} -f {input.genome} --check-ref w --rm-dup exact --atomize)"
        "2> {log}"


rule stratify_by_truth_region:
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


rule bwa_index:
    input:
        "resources/genomes/main.fasta",
    output:
        idx=multiext("resources/genomes/main", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa-index.log",
    # params:
    #     prefix=get_io_prefix(lambda input, output: output[0]),
    wrapper:
        "v1.3.1/bio/bwa/index"


rule bwa_mem:
    input:
        reads=get_fastqs,
        idx=rules.bwa_index.output,
    output:
        temp("results/benchmarking/read-alignments/{sample}.bam"),
    log:
        "logs/bwa-mem/{sample}.log",
    params:
        # index=get_io_prefix(lambda input, output: input.index[0]),
        sorting="samtools",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="coordinate",  # Can be 'queryname' or 'coordinate'.
    threads: 4
    wrapper:
        "v1.3.1/bio/bwa/mem"


rule mark_duplicates:
    input:
        "results/benchmarking/read-alignments/{sample}.bam",
    output:
        bam=temp("results/benchmarking/read-alignments/{sample}.dedup.bam"),
        metrics=temp("results/benchmarking/read-alignments/{sample}.dedup.metrics.txt"),
    log:
        "logs/picard-dedup/{sample}.log",
    params:
        extra="--REMOVE_DUPLICATES true",
    resources:
        mem_mb=1024,
    wrapper:
        "v1.3.1/bio/picard/markduplicates"


rule samtools_index:
    input:
        "results/benchmarking/read-alignments/{sample}.dedup.bam",
    output:
        temp("results/benchmarking/read-alignments/{sample}.dedup.bam.bai"),
    log:
        "logs/samtools-index/{sample}.log",
    wrapper:
        "v1.3.1/bio/samtools/index"


rule mosdepth:
    input:
        bam="results/benchmarking/read-alignments/{sample}.dedup.bam",
        bai="results/benchmarking/read-alignments/{sample}.dedup.bam.bai",
    output:
        temp("results/benchmarking/coverage/{sample}.mosdepth.global.dist.txt"),
        temp("results/benchmarking/coverage/{sample}.quantized.bed.gz"),  # optional, needs to go with params.quantize spec
        summary=temp("results/benchmarking/coverage/{sample}.mosdepth.summary.txt"),  # this named output is required for prefix parsing
    log:
        "logs/mosdepth/{sample}.log",
    params:
        extra="--no-per-base --mapq 59",  # we do not want low MAPQ regions end up being marked as high coverage
        quantize=get_mosdepth_quantize(),
    wrapper:
        "v1.3.1/bio/mosdepth"


rule stratify_by_coverage:
    input:
        variant_call="results/benchmarking/variant-calls/{source}/{workflow}/{sample}.vcf.gz",
        coverage="results/benchmarking/coverage/{sample}.quantized.bed.gz",
    output:
        "results/benchmarking/variant-calls/covered/{source}/{workflow}/{sample}.cov-{cov}.vcf.gz",
    params:
        cov_label=get_cov_label,
    log:
        "logs/stratify-regions/{workflow}/{source}/{sample}/{cov}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "(bedtools intersect -header"
        " -a {input.variant_call}"
        " -b <( zcat {input.coverage} | grep -w  '{params.cov_label}' | bedtools merge -i /dev/stdin ) |"
        " bcftools view -O z -o {output} /dev/stdin)"
        "2>{log}"


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
        truth="results/benchmarking/variant-calls/covered/{source}/{workflow_A}/{sample}.cov-{cov}.vcf.gz",
        truth_idx="results/benchmarking/variant-calls/covered/{source}/{workflow_A}/{sample}.cov-{cov}.vcf.gz.csi",
        query="results/benchmarking/variant-calls/covered/{source}/{workflow_B}/{sample}.cov-{cov}.vcf.gz",
        genome="resources/genomes/main.fasta",
        genome_index="resources/genomes/main.fasta.fai",
    output:
        multiext(
            "results/benchmarking/happy/{source}/{workflow_A}-vs-{workflow_B}/{sample}/cov-{cov}/report",
            ".runinfo.json",
            ".vcf.gz",
            ".summary.csv",
            ".extended.csv",
            ".metrics.json.gz",
            ".roc.all.csv.gz",
        ),
        # ".roc.Locations.SNP.csv.gz",
        # ".roc.Locations.INDEL.csv.gz", # bc of empty vcfs files
        # ".roc.Locations.INDEL.PASS.csv.gz",
    params:
        prefix=lambda w, input, output: output[0].split(".")[0],
        engine="vcfeval",
    log:
        "logs/happy/{source}/{workflow_A}-vs-{workflow_B}/{sample}.cov-{cov}.log",
    wrapper:
        "v1.0.0/bio/hap.py/hap.py"


checkpoint get_samples_with_multiallelic_calls:
    input:
        vcfs=get_benchmark_path(
            "results/benchmarking/variant-calls/covered/{{source}}/{workflow}/{sample}.cov-{{cov}}.vcf.gz",
            remove="sanger",
        ),
        index=get_benchmark_path(
            "results/benchmarking/variant-calls/covered/{{source}}/{workflow}/{sample}.cov-{{cov}}.vcf.gz.csi",
            remove="sanger",
        ),
    output:
        "results/benchmarking/tabels/multiallelic_{source}_{cov}_calls.tsv",
    log:
        "logs/get_samples_with_multiallelic_calls/{source}.cov-{cov}.log",
    conda:
        "../envs/python.yaml"
    params:
        metadata=get_benchmark_path("{workflow},{sample}", remove="sanger"),
    script:
        "../scripts/get_samples_with_multiallelic_calls.py"


checkpoint extract_truth_without_calls:
    input:
        truth=expand(
            "results/benchmarking/variant-calls/covered/{{source}}/sanger/{sample}.cov-{{cov}}.vcf.gz",
            sample=get_samples(),
        ),
        idx=expand(
            "results/benchmarking/variant-calls/covered/{{source}}/sanger/{sample}.cov-{{cov}}.vcf.gz.csi",
            sample=get_samples(),
        ),
    output:
        "results/benchmarking/tabels/truth_without_calls_{source}_{cov}_calls.tsv",
    log:
        "logs/get_truth_with_calls/{source}_{cov}.log",
    conda:
        "../envs/python.yaml"
    params:
        samples=get_samples(),
    script:
        "../scripts/extract_truth_with_calls.py"


rule agg_happy:
    input:
        get_happy_output(
            "results/benchmarking/happy/{{source}}/sanger-vs-{workflow}/{sample}/cov-{{cov}}/report.summary.csv"
        ),
    output:
        "results/benchmarking/tabels/workflow-comparison.source-{source}.cov-{cov}.tsv",
    log:
        "logs/agg_happy/{source}.cov-{cov}.log",
    conda:
        "../envs/python.yaml"
    params:
        metadata=get_happy_output("{workflow},{sample}"),
        platforms=get_happy_platform_data(remove="sanger"),
    script:
        "../scripts/agg_happy.py"


rule agg:
    input:
        expand(
            "results/benchmarking/tabels/workflow-comparison.source-{source}.cov-{cov}.tsv",
            source=["stratified", "normalized-variants"],
            cov=COVERAGES.keys(),
        ),


rule plot_precision_recall:
    input:
        "results/benchmarking/tabels/workflow-comparison.source-{source}.cov-{cov}.tsv",
    output:
        plot="results/benchmarking/plots/variant-calls/precision-recall.source-{source}.cov-{cov}.svg",
        data="results/benchmarking/tabels/variant-calls/precision-recall.source-{source}.cov-{cov}.tsv",
    log:
        "logs/plot_precision_recall.source-{source}.cov-{cov}.log",
    conda:
        "../envs/python.yaml"
    params:
        cov_label=get_cov_label,
    script:
        "../scripts/plot_precision_recall.py"


rule plot_mismatches:
    input:
        "results/benchmarking/tabels/workflow-comparison.stratified.source-{source}.cov-{cov}.tsv",
    output:
        plot="results/benchmarking/plots/variant-missmatches/source-{source}.cov-{cov}.svg",
        data="results/benchmarking/tabels/variant-missmatches/source-{source}.cov-{cov}.tsv",
    log:
        "logs/plot_mismatches.source-{source}.cov-{cov}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_mismatches.py"
