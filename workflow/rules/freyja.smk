rule freyja_variants:
    input:
        bam="results/{date}/read-clipping/hardclipped/pe/{sample}/{sample}.bam",
        ref="resources/genomes/MN908947.fasta",
    output:
        var="results/{date}/lineage-abundance/freyja/{sample}.variants.tsv",
        depth="results/{date}/lineage-abundance/freyja/{sample}.depths.txt",
    log:
        "logs/{date}/freyja/variants/{sample}.log"
    params:
        var=lambda w, output: os.path.splitext(output.var)[0],
    conda:
        "../envs/freyja.yaml"
    shell:
        "freyja variants {input.bam} --variants {params.var} --depths {output.depth} --ref {input.ref} > {log} 2>&1"


rule freyja_demix:
    input:
        var="results/{date}/lineage-abundance/freyja/{sample}.variants.tsv",
        depth="results/{date}/lineage-abundance/freyja/{sample}.depths.txt",
    output:
        "results/{date}/lineage-abundance/freyja/{sample}.demix.txt",
    log:
        "logs/{date}/freyja/demix/{sample}.log"
    conda:
        "../envs/freyja.yaml"
    shell:
        "freyja demix {input.var} {input.depth} --output {output} > {log} 2>&1"


rule aggregate_freyja:
    input:
        demix=lambda wildcards: expand(
            "results/{{date}}/lineage-abundance/freyja/{sample}.demix.txt",
            sample=get_timeseries_samples("sample"),
        ),
    output:
        all="results/{date}/lineage-abundance/freyja/all.demix.csv",
        all_count="results/{date}/lineage-abundance/freyja/all.count.csv",
        pivot="results/{date}/lineage-abundance/freyja/all.pivot.csv",
    log:
        "logs/{date}/aggregate_freyja/all.log",
    params:
        sample=lambda wildcards: get_timeseries_samples("sample"),
        location=lambda wildcards: get_timeseries_samples("location"),
        timestamp=lambda wildcards: get_timeseries_samples("timestamp"),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/aggregate-freyja.py"


rule get_aggregated_freyja:
    input:
        "results/2022-12-21/lineage-abundance/freyja/all.demix.csv",