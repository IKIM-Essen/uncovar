use rule minimap2_bam_sanger as align_sequence with:
    input:
        query="results/benchmarking/backups/{output_type}/{tech}/{workflow}/{sample}",
        target="resources/genomes/main.fasta",
    output:
        "results/benchmarking/sequences/aligned/{output_type}/{tech}/{workflow}/{sample}.bam",
    log:
        "logs/benchmarking/align_sequences/{output_type}/{tech}/{workflow}/{sample}.log",
    params:
        extra="",
        sorting="coordinate",
        sort_extra="",


rule benchmark_sequence:
    input:
        fasta="results/benchmarking/backups/{output_type}/{tech}/{workflow}/{sample}",
        bam="results/benchmarking/sequences/aligned/{output_type}/{tech}/{workflow}/{sample}.bam",
        reference="resources/genomes/main.fasta",
    output:
        "results/benchmarking/sequences/quast/{output_type}/{tech}/{workflow}/{sample}/report.tsv",
    params:
        outdir=get_output_dir,
    log:
        "logs/benchmarking/quast/{output_type}/{tech}/{workflow}/{sample}.log",
    conda:
        "../envs/quast.yaml"
    threads: 2
    shell:
        "quast.py --min-contig 1 --threads {threads} -o {params.outdir} -r {input.reference} "
        "--bam {input.bam} {input.fasta} "
        "> {log} 2>&1"


rule agg_quast:
    input:
        get_output_of_pipelines(
            path="results/benchmarking/backups/{output_type}/{tech}/{workflow}/{sample}.fasta",
            output="consensus",
        ),
    output:
        "results/benchmarking/tabels/quast-sequences.tsv",
    log:
        "logs/agg_quast.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/agg_quast.py"
