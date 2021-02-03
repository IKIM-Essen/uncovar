rule simulate_strain_reads:
    input:
        "resources/genomes/{accession}.fasta",
    output:
        left="resources/benchmarking/{accession}/reads.1.fastq.gz",
        right="resources/benchmarking/{accession}/reads.2.fastq.gz",
    log:
        "logs/mason/{accession}.log",
    conda:
        "../envs/mason.yaml"
    shell:  # median reads in data: 584903
        "mason_simulator -ir {input} -n 584903 -o {output.left} -or {output.right} 2> {log}"


rule test_benchmark_results:
    input:
        get_benchmark_results,
    output:
        "results/benchmarking/strain-calling.csv",
    params:
        true_accessions=get_strain_accessions,
    log:
        "logs/test-benchmark-results.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/test-benchmark-results.py.ipynb"


rule test_assembly_results:
    input:
        "resources/genomes/{accession}.fasta",
        "results/polished-contigs/benchmark-sample-{accession}.fasta",
    output:
        "results/benchmarking/assembly/{accession}.bam",
    log:
        "logs/test-assembly-results/{accession}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 --MD --eqx -ax asm5 {input} -o {output} 2> {log}"


rule summarize_assembly_results:
    input:
        bams=get_assembly_comparisons(bams=True),
        refs=get_assembly_comparisons(bams=False),
    output:
        "results/benchmarking/assembly.csv",
    log:
        "logs/assembly/assembly-results.log",
    conda:
        "../envs/pysam.yaml"
    notebook:
        "../notebooks/assembly-benchmark-results.py.ipynb"
