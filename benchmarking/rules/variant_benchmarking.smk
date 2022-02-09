# TODO check if the pipelines can run on other refernce genome
rule rename_contig:
    input:
        lambda w: get_vcf_of_workflow(w.workflow, w),
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
        "sed 's/{params.search_string}/{params.replace_string}/' {input} > {output} 2> {log}"


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
        truth="results/benchmarking/variant-calls/normalized-variants/{workflow_A}/{sample}.vcf.gz",
        truth_idx="results/benchmarking/variant-calls/normalized-variants/{workflow_A}/{sample}.vcf.gz.csi",
        query="results/benchmarking/variant-calls/normalized-variants/{workflow_B}/{sample}.vcf.gz",
        genome="resources/genomes/main.fasta",
        genome_index="resources/genomes/main.fasta.fai",
    output:
        multiext(
            "results/benchmarking/happy/{workflow_A}-vs-{workflow_B}/{sample}/report",
            ".runinfo.json",
            ".vcf.gz",
            ".summary.csv",
            ".extended.csv",
            ".metrics.json.gz",
            ".roc.all.csv.gz",
            ".roc.Locations.INDEL.csv.gz",
            ".roc.Locations.INDEL.PASS.csv.gz",
            ".roc.Locations.SNP.csv.gz",
        ),
    params:
        prefix=lambda w, input, output: output[0].split(".")[0],
        engine="vcfeval",
    log:
        "logs/happy/{workflow_A}-vs-{workflow_B}/{sample}.log",
    wrapper:
        "v1.0.0/bio/hap.py/hap.py"


rule agg_normalize_calls:
    input:
        lambda w: expand(
            "results/benchmarking/variant-calls/normalized-variants/{workflow}/{sample}.vcf.gz",
            workflow=PIPELINES["nanopore"],
            sample=get_nanopore_samples(w),
        ),
        lambda w: expand(
            "results/benchmarking/variant-calls/normalized-variants/{workflow}/{sample}.vcf.gz",
            workflow=PIPELINES["illumina"],
            sample=get_illumina_samples(w),
        ),


rule agg_happy:
    input:
        lambda w: expand(
            "results/benchmarking/happy/sanger-vs-{workflow}/{sample}/report.runinfo.json",
            workflow=PIPELINES["nanopore"],
            sample=get_nanopore_samples(w),
        ),
        lambda w: expand(
            "results/benchmarking/happy/sanger-vs-{workflow}/{sample}/report.runinfo.json",
            workflow=PIPELINES["illumina"],
            sample=get_illumina_samples(w)[0],
        ),


# rule benchmark_variants:
#     input:
#         truth=lambda w: get_vcf_of_workflow(w.pipeline_1, w),  # sanger vcf
#         query=lambda w: get_vcf_of_workflow(w.pipeline_2, w),  # variant calls
#         truth_regions="results/benchmarking/sanger/aligned/{sample}.bed",  # sanger bed
#         genome="resources/genomes/main.fasta",
#         genome_index="resources/genomes/main.fasta.fai",
#     output:
#         multiext(
#             "results/benchmarking/happy/happy/truth~{pipeline_1}-vs-{pipeline_2}/{sample}/{sample}",
#             ".runinfo.json",
#             ".vcf.gz",
#             ".summary.csv",
#             ".extended.csv",
#             ".metrics.json.gz",
#             ".roc.all.csv.gz",
#             ".roc.Locations.INDEL.csv.gz",
#             ".roc.Locations.INDEL.PASS.csv.gz",
#             ".roc.Locations.SNP.csv.gz",
#             ".roc.tsv",
#         ),
#     params:
#         engine="vcfeval",
#         prefix=lambda wc, input, output: output[0].split(".")[0],
#         ## parameters such as -L to left-align variants
#         extra="--verbose",
#     log:
#         "logs/happy/truth~{pipeline_1}-vs-{pipeline_2}/{sample}.log",
#     threads: 2
#     wrapper:
#         "v1.0.0/bio/hap.py/hap.py"
# rule agg_happy:
#     input:
#         lambda w: expand(
#             "results/benchmarking/happy/happy/truth~sanger-vs-{pipeline}/{sample}/{sample}.runinfo.json",
#             pipeline=PIPELINES["nanopore"],
#             sample=get_nanopore_samples(w),
#         ),
#         lambda w: expand(
#             "results/benchmarking/happy/happy/truth~sanger-vs-{pipeline}/{sample}/{sample}.runinfo.json",
#             pipeline=PIPELINES["illumina"],
#             sample=get_illumina_samples(w),
#         ),
