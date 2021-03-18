rule freebayes:
    input:
        ref=get_reference(),
        ref_idx=get_reference(".fai"),
        # you can have a list of samples here
        samples="results/{date}/recal/ref~{reference}/{sample}.bam",
        index="results/{date}/recal/ref~{reference}/{sample}.bam.bai",
    output:
        temp("results/{date}/candidate-calls/ref~{reference}/{sample}.bcf"),
    log:
        "logs/{date}/freebayes/ref~{reference}/{sample}.log",
    params:
        # genotyping is performed by varlociraptor, hence we deactivate it in freebayes by 
        # always setting --pooled-continuous
        extra=(
            "--pooled-continuous --min-alternate-count 1 --min-alternate-fraction 0.01"
        ),
    wrapper:
        "0.68.0/bio/freebayes"


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


rule varlociraptor_preprocess:
    input:
        ref=get_reference(),
        ref_idx=get_reference(".fai"),
        candidates="results/{date}/candidate-calls/ref~{reference}/{sample}.bcf",
        bam="results/{date}/recal/ref~{reference}/{sample}.bam",
        bai="results/{date}/recal/ref~{reference}/{sample}.bam.bai",
    output:
        "results/{date}/observations/ref~{reference}/{sample}.bcf",
    params:
        depth=config["variant-calling"]["max-read-depth"],
    log:
        "logs/{date}/varlociraptor/preprocess/ref~{reference}/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --max-depth {params.depth} --output {output} 2> {log}"


rule varlociraptor_call:
    input:
        obs="results/{date}/observations/ref~{reference}/{sample}.bcf",
        scenario="results/{date}/scenarios/{sample}.yaml",
    output:
        temp("results/{date}/calls/ref~{reference}/{sample}.bcf"),
    params:
        biases=get_varlociraptor_bias_flags,
    log:
        "logs/{date}/varlociraptor/call/ref~{reference}/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants {params.biases} generic --obs {wildcards.sample}={input.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"
