# call candidate variants with random index
# alternative samtools mpileup / bcftools mpileup
# freebayes looks at short haplotypes, able to detect MNV -> more informative
rule freebayes:
    input:
        ref=get_reference(),
        ref_idx=get_reference(".fai"),
        # you can have a list of samples here
        samples="results/recal/ref~{reference}/{sample}.bam",
        index="results/recal/ref~{reference}/{sample}.bam.bai",
    output:
        "results/candidate-calls/ref~{reference}/{sample}.bcf",
    log:
        "logs/freebayes/ref~{reference}/{sample}.log",
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
        report(
            "results/scenarios/{sample}.yaml",
            caption="../report/scenario.rst",
            category="Variant calling scenarios",
        ),
    log:
        "logs/render-scenario/{sample}.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "sed 's/sample:/{wildcards.sample}:/' {input} > {output}"


rule varlociraptor_preprocess:
    input:
        ref=get_reference(),
        ref_idx=get_reference(".fai"),
        candidates="results/candidate-calls/ref~{reference}/{sample}.bcf",
        bam="results/recal/ref~{reference}/{sample}.bam",
        bai="results/recal/ref~{reference}/{sample}.bam.bai",
    output:
        "results/observations/ref~{reference}/{sample}.bcf",
    params:
        depth=config["variant-calling"]["max-read-depth"],
    log:
        "logs/varlociraptor/preprocess/ref~{reference}/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants --candidates {input.candidates} "
        "{input.ref} --bam {input.bam} --max-depth {params.depth} --output {output} 2> {log}"


rule varlociraptor_call:
    input:
        obs="results/observations/ref~{reference}/{sample}.bcf",
        scenario="results/scenarios/{sample}.yaml",
    output:
        "results/calls/ref~{reference}/{sample}.bcf",
    log:
        "logs/varlociraptor/call/ref~{reference}/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor "
        "call variants generic --obs {wildcards.sample}={input.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"
