rule assembly:
    input:
        fastq1="results/trimmed/{sample}.1.fastq.gz",
        fastq2="results/trimmed/{sample}.2.fastq.gz",
    output:
        "results/assembly/{sample}/final.contigs.fa",
    log:
        "logs/megahit/{sample}.log",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/megahit.yaml"
    shell:
        "rm -r {params.outdir}; megahit -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir} 2> {log}"


rule order_contigs:
    input:
        contigs="results/assembly/{sample}/final.contigs.fa",
        reference="resources/genomes/main.fasta",
    output:
        "results/ordered_contigs/{sample}/ragoo_output/ragoo.fasta",
    log:
        "logs/ragoo/{sample}.log",
    params:
        outdir=lambda x, output: os.path.dirname(os.path.dirname(output[0])),
    threads: 8
    conda:
        "../envs/ragoo.yaml"
    shell:  # currently there is no conda package for mac available. Manuell download via https://github.com/malonge/RaGOO
        "(mkdir -p {params.outdir} && cd {params.outdir} && "
        "ragoo.py ../../../{input.contigs} ../../../{input.reference}) 2> {log}"

# rule align_ordered_contigs:
#     input:
#         "resources/genomes/main.fasta",
#         "results/ordered_contigs/{sample}/ragoo_output/ragoo.fasta",
#     output:
#         "results/ordered_contigs/{sample}/{sample}.bam",
#     log:
#         "logs/test_ragoo/{sample}.log",
#     conda:
#         "../envs/minimap2.yaml"
#     shell:
#         "minimap2 -ax asm5 {input} -o {output} 2> {log}"

# Which contigs to vis? ordered or unordered?
# rule quast:
#     input:
#         reference="resources/genomes/main.fasta",
#         bams=get_all_ordered_samples_as_bams,
#         fastas=get_all_ordered_samples_as_fastas
#     output:
#         directory('results/quast')
#     log:
#         "logs/quast.log"
#     params:
#         bam_list = lambda x, input: ",".join(input.bams)
#     conda:
#         "../envs/quast.yaml"
#     threads: 8
#     shell:
#         "(quast.py --threads {threads} -o {output} -r {input.reference} --eukaryote --bam {params.bam_list} {input.fastas}) 2> {log}"


# TODO add plot that visualizes assembly quality
# TODO blast smaller contigs to determine contamination?
