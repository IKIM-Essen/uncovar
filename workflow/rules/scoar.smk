# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


configfile: "config/config.yaml"


rule collect_minority_mutations:
    input:
        bcfs=lambda wildcards: expand(
            "results/{{date}}/annotated-calls/ref~main/annot~orf/{sample}.bcf",
            sample=get_samples_for_date(wildcards),
        ),
        pan_calls="results/{date}/tables/pangolin_calls_per_stage.csv",
    output:
        "results/{date}/scoar/collection_minority_variants.csv",
    log:
        "logs/{date}/scoar/collectminor.log",
    conda:
        "../envs/collectminor.yaml"
    script:
        "../scripts/collect_minority_mutations.py"


rule scoar:
    input:
        "results/{date}/scoar/collection_minority_variants.csv",
    output:
        "results/{date}/scoar/scoar_resistogram.csv",
    log:
        "logs/{date}/scoar/scoar.log",
    conda:
        "../envs/scoar.yaml"
    shell:
        "../../resources/scoar.py {input} -p -o {output} 2> {log}"



