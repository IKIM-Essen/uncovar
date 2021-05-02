rule count_mapped_reads:
    input:
        bam="results/{date}/mapped/ref~main/{sample}.bam",
        fasta="resources/genomes/main.fasta",
    output:
        "results/{date}/savage/ref~main/{sample}.split",
    log:
        "logs/{date}/count-mapped/{sample}.log",
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -c -F 260 {input} > {output} 2> {log}"


rule savage:
    input:
        ref="resources/genomes/main.fasta",
        idx="resources/genomes/main.fasta.fai",
        fastq1="results/{date}/clipped-reads/{sample}.1.fastq.gz",
        fastq2="results/{date}/clipped-reads/{sample}.2.fastq.gz",
        bam="results/{date}/mapped/ref~main/{sample}.bam",
    output:
        assembly=directory("results/{date}/savage/{sample}"),
        depth="results/{date}/depth/{sample}.txt",
    params:
        denovo=True,
    conda:
        "../envs/savage.yaml"
    threads: 16
    script:
        "../scripts/savage.py"


rule deploy_virus_vg:
    output:
        directory("resources/virus-vg"),
    log:
        "logs/deploy-virus-vg.log",
    conda:
        "../envs/git.yaml"
    shell:
        "(git clone https://bitbucket.org/jbaaijens/virus-vg.git {output}; "
        "cd {output}; git checkout 6eec5149711cff0682f737b348c3dde646973cd1) 2> {log}"


rule build_variation_graph:
    input:
        fastq1="results/{date}/clipped-reads/{sample}.1.fastq.gz",
        fastq2="results/{date}/clipped-reads/{sample}.2.fastq.gz",
        savage="results/{date}/savage/{sample}",
        virus_vg="resources/virus-vg",
    output:
        "results/virus-vg/variation-graph/{sample}/contig_graph.final.gfa",
        "results/virus-vg/variation-graph/{sample}/node_abundances.txt",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    log:
        "logs/{date}/virus-vg/build-variation-graph/{sample}.log",
    conda:
        "../envs/virus-vg.yaml"
    threads: 16
    shell:
        "cd {params.outdir}; {input.virus_vg}/scripts/build_graph_msga.py -f {input.fastq1} -r {input.fastq2} "
        "-c {input.savage}/contigs_stage_c.fasta -vg $(which vg) -t {threads} 2> {log}"


rule build_haplotypes:
    input:
        gfa="results/virus-vg/variation-graph/{sample}/contig_graph.final.gfa",
        abundances="results/virus-vg/variation-graph/{sample}/node_abundances.txt",
        virus_vg="resources/virus-vg",
        depth="results/{date}/depth/{sample}.txt",
    output:
        "results/{date}/assembly/virus-vg/{sample}.contigs.fasta",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    script:
        "../scripts/virus-vg-optimize-strains.py"
