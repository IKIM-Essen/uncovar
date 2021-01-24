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
    shell:
        "mason_simulator -ir {input} -n 30424 -o {output.left} -or {output.right} 2> {log}"


rule test_benchmark_results:
    input:
        get_benchmark_results,
    output:
        "results/benchmarking.csv",
    params:
        true_accessions=get_strain_accessions,
    log:
        "logs/test-benchmark-results.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/test-benchmark-results.py.ipynb"
