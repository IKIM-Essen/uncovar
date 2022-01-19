# source: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
rule artic_guppyplex:
    input:
        get_barcode_path,
    output:
        "results/benchmarking/other-pipelines/artic/guppyplex/{sample}.fasta",
    log:
        "logs/artic_guppyplex/{sample}.log",
    conda:
        "../../envs/artic.yaml"
    shell:
        "artic guppyplex --min-length 400 --max-length 700 "
        "--directory {input} --output {output} > {log}"


# rule artic_minion:
#     input:
#         "results/benchmarking/other-pipelines/artic/guppyplex/{sample}.fasta",
#     output:
#         "",
#     log:
#         "logs/artic_minion.log"
#     conda:
#         "../../envs/.yaml"
#     shell:
#         ""
