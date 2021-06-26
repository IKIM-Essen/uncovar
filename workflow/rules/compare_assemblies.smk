rule assembly_megahit:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/megahit/{sample}/{sample}.contigs.fasta",
    log:
        "logs/{date}/megahit/{sample}.log",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/megahit.yaml"
    shell:
        "(megahit -1 {input.fastq1} -2 {input.fastq2} --out-dir {params.outdir} -f && "
        "mv {params.outdir}/final.contigs.fa {output.contigs} ) > {log} 2>&1"


rule assembly_trinity:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        "results/{date}/assembly/trinity/{sample}/{sample}.contigs.fasta",
    log:
        'logs/{date}/trinity/{sample}.log'
    params:
        extra=""
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads: 8
    shell:
        "(Trinity --left {input.fastq1} --right {input.fastq2} --CPU {snakemake.threads} --seqType fq --output {params.outdir} && "
        "mv {params.outdir}/Trinity.fasta {output.contigs} ) > {log} 2>&1"


rule assembly_metaspades:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/metaspades/{sample}/{sample}.contigs.fasta",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    log:
        "logs/{date}/metaSPAdes/{sample}.log",
    threads: 8
    conda:
        "../envs/spades.yaml"
    shell:
        "(metaspades.py -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir} -t {threads} && "
        "mv {params.outdir}/contigs.fasta {output.contigs}) > {log} 2>&1"


rule assembly_coronaspades:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/coronaspades/{sample}/{sample}.contigs.fasta",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    log:
        "logs/{date}/coronaSPAdes/{sample}.log",
    threads: 8
    conda:
        "../envs/spades.yaml"
    shell:
        "(coronaspades.py -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir} -t {threads} && "
        "mv {params.outdir}/contigs.fasta {output.contigs}) > {log} 2>&1"


rule assembly_spades:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/spades/{sample}/{sample}.contigs.fasta",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    log:
        "logs/{date}/SPAdes/{sample}.log",
    threads: 8
    conda:
        "../envs/spades.yaml"
    shell:
        "(spades.py -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir} -t {threads} && "
        "mv {params.outdir}/contigs.fasta {output.contigs}) > {log} 2>&1"


rule assembly_rnaviralspades:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/rnaviralspades/{sample}/{sample}.contigs.fasta",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    log:
        "logs/{date}/rnaviralSPAdes/{sample}.log",
    threads: 8
    conda:
        "../envs/spades.yaml"
    shell:
        "(rnaviralspades.py -1 {input.fastq1} -2 {input.fastq2} -o {params.outdir} -t {threads} && "
        "mv {params.outdir}/contigs.fasta {output.contigs}) > {log} 2>&1"


rule order_contigs:
    input:
        contigs=get_contigs,
        reference="resources/genomes/main.fasta",
    output:
        temp("results/{date}/contigs/ordered-unfiltered/{sample}.fasta"),
    log:
        "logs/{date}/ragoo/{sample}.log",
    params:
        outdir=lambda x, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/ragoo.yaml"
    shell:  # currently there is no conda package for mac available. Manuell download via https://github.com/malonge/RaGOO
        'if [ -d "{params.outdir}/{wildcards.sample}" ]; then rm -Rf {params.outdir}/{wildcards.sample}; fi && '
        "(mkdir -p {params.outdir}/{wildcards.sample} && cd {params.outdir}/{wildcards.sample} && "
        "ragoo.py ../../../../../{input.contigs} ../../../../../{input.reference} && "
        "cd ../../../../../ && mv {params.outdir}/{wildcards.sample}/ragoo_output/ragoo.fasta {output}) > {log} 2>&1"


rule filter_chr0:
    input:
        "results/{date}/contigs/ordered-unfiltered/{sample}.fasta",
    output:
        "results/{date}/contigs/ordered/{sample}.fasta",
    log:
        "logs/{date}/ragoo/{sample}_cleaned.log",
    params:
        sample=lambda wildcards: wildcards.sample,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/ragoo-remove-chr0.py"