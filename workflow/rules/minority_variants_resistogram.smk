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
        "results/{date}/resistogram/collection_minority_variants.csv",
    log:
        "logs/{date}/resistogram/collectminor.log",
    conda:
        "../envs/collectminor.yaml"
    script:
        "../scripts/collect_minority_mutations.py"


rule generate_resistogram:
    input:
        escaping_mutations="resources/resistogram/mabs.json",
        factors="resources/resistogram/factortab.csv",
        allmutationsfound="results/{date}/resistogram/collection_minority_variants.csv",
    output:
        "results/{date}/resistogram/resistogram.json",
    log:
        "logs/{date}/resistogram/resistogram.log",
    conda:
        "../envs/resistogram.yaml"
    script:
        "../scripts/generate_resistogram.py"


# rule post_to_web_ui:
#     input:
#         "results/{{date}}/resistogram/resistogram.json"
#     output:
#         temp("results/{{date}}/resistogram/successfulpost.txt")
#     log:
#         "logs/{{date}}/resistogram/post.log"
#     run:
#         import requests
#         url = 'https://www.fancyUI.com/covid/demopage.php'
#         obj = {'resistogram': 'results/{{date}}/resistogram/resistogram.json'}
#         requests.post(url, json = obj)
