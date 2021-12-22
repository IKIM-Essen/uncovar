# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule freebayes:
    input:
        ref=get_reference(),
        ref_idx=get_reference(".fai"),
        # you can have a list of samples here
        samples="results/{date}/recal/ref~{reference}/{sample}.bam",
        index="results/{date}/recal/ref~{reference}/{sample}.bam.bai",
    output:
        temp("results/{date}/candidate-calls/ref~{reference}/{sample}.small.bcf"),
    params:
        # genotyping is performed by varlociraptor, hence we deactivate it in freebayes by
        # always setting --pooled-continuous
        extra=(
            "--pooled-continuous --min-alternate-count 1 --min-alternate-fraction 0.01"
        ),
    log:
        "logs/{date}/freebayes/ref~{reference}/{sample}.log",
    wrapper:
        "0.80.1/bio/freebayes"


# TODO check delly single end mode
rule delly:
    input:
        ref=get_reference(),
        ref_idx=get_reference(".fai"),
        sample="results/{date}/recal/ref~{reference}/{sample}.bam",
        sample_idx="results/{date}/recal/ref~{reference}/{sample}.bam.bai",
    output:
        temp("results/{date}/candidate-calls/ref~{reference}/{sample}.structural.bcf"),
    log:
        "logs/{date}/delly/ref~{reference}/{sample}.log",
    conda:
        "../envs/delly.yaml"
    script:
        "../scripts/delly.py"


rule medaka_variant:
    input:
        ref=get_reference(),
        sample="results/{date}/recal/ref~{reference}/{sample}.bam",
        bai="results/{date}/recal/ref~{reference}/{sample}.bam.bai",
    output:
        temp(
            "results/{date}/candidate-calls/ref~{reference}/{sample}.homopolymer-medaka.vcf"
        ),
    params:
        outdir=get_output_dir,
    log:
        "logs/{date}/medaka/variant/ref~{reference}/{sample}.log",
    conda:
        "../envs/medaka.yaml"
    threads: 4
    shell:
        "(medaka_variant -i {input.sample} -f {input.ref} -o {params.outdir}/{wildcards.sample} -t {threads} &&"
        " mv $(cat {log}| grep '\- Final VCF written to' | sed s/'- Final VCF written to '//\ ) {output})"
        " > {log} 2>&1"


rule longshot:
    input:
        ref=get_reference(),
        bam="results/{date}/recal/ref~{reference}/{sample}.bam",
        bai="results/{date}/recal/ref~{reference}/{sample}.bam.bai",
    output:
        temp(
            "results/{date}/candidate-calls/ref~{reference}/{sample}.homopolymer-longshot.vcf"
        ),
    params:
        reference_name=lambda w: config["virus-reference-genome"]
        if w.reference == "main"
        else f"{w.reference}.1",
    log:
        "logs/{date}/longshot/ref~{reference}/{sample}.log",
    conda:
        "../envs/longshot.yaml"
    shell:
        "(longshot -P 0 -F -A --no_haps --bam {input.bam} --ref {input.ref} --out {output} &&"
        " sed -i '2 i\##contig=<ID={params.reference_name}>' {output})"
        " 2> {log}"


# TODO adjust for ion


rule vcf_2_bcf:
    input:
        "results/{date}/candidate-calls/ref~{reference}/{sample}.{varrange}.vcf",
    output:
        temp("results/{date}/candidate-calls/ref~{reference}/{sample}.{varrange}.bcf"),
    log:
        "logs/{date}/vcf_2_bcf/ref~{reference}/{sample}.{varrange}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools view -O u -o {output} {input} 2> {log}"


rule render_scenario:
    input:
        local(get_resource("scenario.yaml")),
    output:
        "results/{date}/scenarios/{sample}.yaml",
    log:
        "logs/{date}/render-scenario/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "sed 's/sample:/{wildcards.sample}:/' {input} > {output}"


rule varlociraptor_alignment_properties:
    input:
        ref=get_reference(),
        ref_idx=get_reference(".fai"),
        bam="results/{date}/recal/ref~{reference}/{sample}.bam",
    output:
        temp("results/{date}/alignment-properties/ref~{reference}/{sample}.json"),
    log:
        "logs/{date}/varlociraptor/estimate-alignment-properties/ref~{reference}/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor estimate alignment-properties {input.ref} --bam {input.bam} > {output} 2> {log}"


rule varlociraptor_preprocess:
    input:
        ref=get_reference(),
        ref_idx=get_reference(".fai"),
        candidates=get_candidate_variants,
        bam="results/{date}/recal/ref~{reference}/{sample}.bam",
        bai="results/{date}/recal/ref~{reference}/{sample}.bam.bai",
    output:
        temp("results/{date}/observations/ref~{reference}/{sample}.{varrange}.bcf"),
    params:
        depth=config["variant-calling"]["max-read-depth"],
    log:
        "logs/{date}/varlociraptor/preprocess/ref~{reference}/{sample}.{varrange}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --max-depth {params.depth} --output {output} 2> {log}"


rule varlociraptor_call:
    input:
        obs="results/{date}/observations/ref~{reference}/{sample}.{varrange}.bcf",
        scenario="results/{date}/scenarios/{sample}.yaml",
    output:
        temp("results/{date}/calls/ref~{reference}/{sample}.{varrange}.bcf"),
    params:
        biases=get_varlociraptor_bias_flags,
    log:
        "logs/{date}/varlociraptor/call/ref~{reference}/{sample}.{varrange}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants {params.biases} generic --obs {wildcards.sample}={input.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"


rule merge_varranges:
    input:
        calls=lambda wildcards: expand(
            "results/{{date}}/calls/ref~{{reference}}/{{sample}}.{varrange}.bcf",
            varrange=get_varrange(wildcards),
        ),
        idx=lambda wildcards: expand(
            "results/{{date}}/calls/ref~{{reference}}/{{sample}}.{varrange}.bcf.csi",
            varrange=get_varrange(wildcards),
        ),
    output:
        "results/{date}/calls/ref~{reference}/{sample}.bcf",
    log:
        "logs/{date}/merge-calls/ref~{reference}/{sample}.log",
    params:
        "-a -Ob",
    wrapper:
        "0.69.0/bio/bcftools/concat"
