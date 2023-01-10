rule collect_minority_mutations: 
    input: 
        bcfs="results/{date}/annotated-calls/ref~main/annot~orf/{sample}.bcf",
        pan_calls="results/{date}/tables/pangolin_calls_per_stage.csv",
    output: 
        "results/{date}/tables/collection_minority_variants.csv",
    log:
        "logs/{date}/resistogram/collectminor.log"
    conda:
        "../envs/collectminor.yaml"
    script:
        "../scripts/collect_minority_mutations.py"

rule generate_resistogram:
    input: 
        escaping_mutations="resources/resistogram/mabs.json",
        factors="resources/resistogram/factortab.csv",
        allmutationsfound="results/{date}/tables/collection_minority_variants.csv",
    output:
        "results/{date}/tables/resistogram.yaml"
    log:
        "logs/{date}/resistogram/resistogram.log"
    conda:
        "../envs/resistogram.yaml"
    script:
        "../scripts/generate_resistogram.py"