rule simulate_strain_reads:
    input:
        get_genome_fasta,
    output:
        left=temp("resources/benchmarking/{accession}/reads.1.fastq.gz"),
        right=temp("resources/benchmarking/{accession}/reads.2.fastq.gz"),
    params:
        no_reads=lambda wildcards: no_reads(wildcards),
    log:
        "logs/mason/benchmarking/{accession}.log",
    conda:
        "../envs/mason.yaml"
    shell:  # median reads in data: 584903
        "mason_simulator -ir {input} -n {params.no_reads} -o {output.left} -or {output.right} 2> {log}"


rule mix_strain_reads:
    input:
        left=expand(
            "resources/benchmarking/{mix}/reads.1.fastq.gz",
            mix=[
                "{MIXTURE_PART_INDICATOR}{{strain_{no}}}".format(
                    MIXTURE_PART_INDICATOR=MIXTURE_PART_INDICATOR, no=i
                )
                for i in range(config["mixtures"]["no_strains"])
            ],
        ),
        right=expand(
            "resources/benchmarking/{mix}/reads.2.fastq.gz",
            mix=[
                "{MIXTURE_PART_INDICATOR}{{strain_{no}}}".format(
                    MIXTURE_PART_INDICATOR=MIXTURE_PART_INDICATOR, no=i
                )
                for i in range(config["mixtures"]["no_strains"])
            ],
        ),
    output:
        left=temp(
            expand(
                "resources/mixtures/{mix}/reads.1.fastq.gz",
                mix="".join(
                    [
                        "{MIXTURE_PART_INDICATOR}{{strain_{no}}}".format(
                            MIXTURE_PART_INDICATOR=MIXTURE_PART_INDICATOR,
                            no=i,
                        )
                        for i in range(config["mixtures"]["no_strains"])
                    ]
                ),
            )
        ),
        right=temp(
            expand(
                "resources/mixtures/{mix}/reads.2.fastq.gz",
                mix="".join(
                    [
                        "{MIXTURE_PART_INDICATOR}{{strain_{no}}}".format(
                            MIXTURE_PART_INDICATOR=MIXTURE_PART_INDICATOR,
                            no=i,
                        )
                        for i in range(config["mixtures"]["no_strains"])
                    ]
                ),
            )
        ),
    log:
        "logs/mix_strain_reads/{}".format(
            "".join(
                [
                    "{MIXTURE_PART_INDICATOR}{{strain_{no}}}".format(
                        MIXTURE_PART_INDICATOR=MIXTURE_PART_INDICATOR, no=i
                    )
                    for i in range(config["mixtures"]["no_strains"])
                ]
            )
        ),
    shell:
        "(zcat {input.left} > {output.left} &&"
        "zcat {input.right} > {output.right}) 2>{log}"


rule test_benchmark_results:
    input:
        get_benchmark_results,
    output:
        "results/benchmarking/strain-calling.csv",
    params:
        true_accessions=get_strain_accessions,
    log:
        "logs/test-benchmark-results.log",
    conda:
        "../envs/python.yaml"
    notebook:
        "../notebooks/test-benchmark-results.py.ipynb"


rule test_assembly_results:
    input:
        "resources/genomes/{accession}.fasta",
        get_assembly_result,
    output:
        "results/benchmarking/assembly/{assembly_type}/{accession}.bam",
    log:
        "logs/test-assembly-results/{assembly_type}/{accession}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 --MD --eqx -ax asm5 {input} -o {output} 2> {log}"


rule summarize_assembly_results:
    input:
        bams=get_assembly_comparisons(bams=True),
        refs=get_assembly_comparisons(bams=False),
    output:
        "results/benchmarking/assembly/{assembly_type}.csv",
    log:
        "logs/summarize-assembly-results/{assembly_type}/assembly-results.log",
    conda:
        "../envs/pysam.yaml"
    notebook:
        "../notebooks/assembly-benchmark-results.py.ipynb"


rule test_non_cov2:
    input:
        pangolin=get_non_cov2_calls(from_caller="pangolin"),
        kallisto=get_non_cov2_calls(from_caller="kallisto"),
    output:
        "results/benchmarking/non-sars-cov-2.csv",
    params:
        accessions=get_non_cov2_accessions(),
    log:
        "logs/benchmarking/summarize_non_cov2.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize-non-cov2.py"


rule report_non_cov2:
    input:
        summary="results/benchmarking/non-sars-cov-2.csv",
        call_plots=expand(
            "results/benchmarking/plots/strain-calls/non-cov2-{accession}.strains.{caller}.svg",
            accession=get_non_cov2_accessions(),
            caller=["pangolin", "kallisto"],
        ),
    output:
        report(
            directory("results/benchmarking/html"),
            htmlindex="index.html",
            category="Test results",
        ),
    log:
        "logs/report_non_cov2.log",
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt csv-report -s '\t' {input.summary} {output}"


checkpoint generate_mixtures:
    input:
        "results/benchmarking/tables/strain-genomes.txt",
    output:
        "results/benchmarking/tables/mixtures.txt",
    params:
        mixtures=generate_mixtures,
    log:
        "logs/generate_mixtures.log",
    script:
        "../scripts/generate-mixtures.py"


rule evaluate_strain_call_error:
    input:
        get_mixture_results,
    output:
        "results/benchmarking/tables/{caller}-strain-call-error.csv",
    params:
        max_reads=config["mixtures"]["max_reads"],
        prefix=MIXTURE_PREFIX,
        separator=MIXTURE_PART_INDICATOR,
        percentage=MIXTURE_PERCENTAGE_INDICATOR,
    log:
        "logs/evaluate-{caller}-strain-call-error.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/evaluate-strain-call-error.py"


rule plot_strain_call_error:
    input:
        "results/benchmarking/tables/{caller}-strain-call-error.csv",
    output:
        "results/benchmarking/plots/{caller}-strain-call-error-heatmap.svg",
        "results/benchmarking/plots/{caller}-strain-call-error-false-predictions.svg",
        "results/benchmarking/plots/{caller}-strain-call-error-content-false-predictions.svg",
    log:
        "logs/plot-{caller}-strain-call-error.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-caller-error.py"


rule assembly_comparison_megahit_std:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/{sample}/megahit-std/{sample}.contigs.fasta",
    log:
        "logs/{date}/megahit-std/{sample}.log",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/megahit.yaml"
    shell:
        "(megahit -1 {input.fastq1} -2 {input.fastq2} --out-dir {params.outdir} -f && "
        "mv {params.outdir}/final.contigs.fa {output.contigs} ) > {log} 2>&1"


rule assembly_comparison_megahit_meta_large:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/{sample}/megahit-meta-large/{sample}.contigs.fasta",
    log:
        "logs/{date}/megahit-meta-large/{sample}.log",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/megahit.yaml"
    shell:
        "(megahit -1 {input.fastq1} -2 {input.fastq2} --presets meta-large --out-dir {params.outdir} -f && "
        "mv {params.outdir}/final.contigs.fa {output.contigs} ) > {log} 2>&1"


rule assembly_comparison_megahit_meta_sensitive:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/{sample}/megahit-meta-sensitive/{sample}.contigs.fasta",
    log:
        "logs/{date}/megahit-meta-sensitive/{sample}.log",
    params:
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/megahit.yaml"
    shell:
        "(megahit -1 {input.fastq1} -2 {input.fastq2} --presets meta-sensitive --out-dir {params.outdir} -f && "
        "mv {params.outdir}/final.contigs.fa {output.contigs} ) > {log} 2>&1"


rule assembly_comparison_trinity:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        "results/{date}/assembly/{sample}/trinity/{sample}.contigs.fasta",
    log:
        "logs/{date}/trinity/{sample}.log",
    params:
        extra="",
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/trinity.yaml"
    shell:
        "(Trinity --left {input.fastq1} --max_memory 16G --right {input.fastq2} --CPU {threads} --seqType fq --output {params.outdir} && "
        "mv {params.outdir}/Trinity.fasta {output} ) > {log} 2>&1"


# rule assembly_comparison_abyss:
#     input:
#         fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
#         fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
#     output:
#         "results/{date}/assembly/{sample}/abyss/{sample}.contigs.fasta",
#     log:
#         "logs/{date}/abyss/{sample}.log",
#     params:
#         extra="",
#         outdir=lambda w, output: os.path.dirname(output[0]),
#     threads: 8
#     conda:
#         "../envs/abyss.yaml"
#     shell:
#         """
#         if [ -d "{params.outdir}" ]; then rm -Rf {params.outdir}; fi
#         (mkdir -p {params.outdir} && cd {params.outdir}
#         abyss-pe np={threads} name={wildcards.sample} k=96 in='../../../../../{input.fastq1} ../../../../../{input.fastq2}'
#         cd ../../../../../ && mv {params.outdir}/{wildcards.sample}-contigs.fa {output} ) > {log} 2>&1
#         """


rule assembly_comparison_velvet:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        "results/{date}/assembly/{sample}/velvet/{sample}.contigs.fasta",
    log:
        "logs/{date}/velvet/{sample}.log",
    params:
        extra="",
        outdir=lambda w, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/velvet.yaml"
    shell:
        """
        velveth {params.outdir} 21 -fastq.gz -shortPaired {input.fastq1} {input.fastq2} > {log} 2>&1
        velvetg {params.outdir} -ins_length 150 -exp_cov 10 >> {log} 2>&1
        mv {params.outdir}/contigs.fa {output} >> {log} 2>&1
        """


rule assembly_comparison_metaspades:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/{sample}/metaspades/{sample}.contigs.fasta",
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


rule assembly_comparison_coronaspades:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/{sample}/coronaspades/{sample}.contigs.fasta",
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


rule assembly_comparison_spades:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs="results/{date}/assembly/{sample}/spades/{sample}.contigs.fasta",
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


rule assembly_comparison_rnaviralspades:
    input:
        fastq1=lambda wildcards: get_reads_after_qc(wildcards, read="1"),
        fastq2=lambda wildcards: get_reads_after_qc(wildcards, read="2"),
    output:
        contigs=(
            "results/{date}/assembly/{sample}/rnaviralspades/{sample}.contigs.fasta"
        ),
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


rule order_contigs_assembly_comparison:
    input:
        contigs="results/{date}/assembly/{sample}/{assembler}/{sample}.contigs.fasta",
        reference="resources/genomes/main.fasta",
    output:
        "results/{date}/assembly/{sample}/{assembler}/{sample}.ordered.contigs.fasta",
    log:
        "logs/{date}/ragoo/{assembler}_{sample}.log",
    params:
        outdir=lambda x, output: os.path.dirname(output[0]),
    threads: 8
    conda:
        "../envs/ragoo.yaml"
    shell:  # currently there is no conda package for mac available. Manuell download via https://github.com/malonge/RaGOO
        "(cd {params.outdir} && "
        "ragoo.py ../../../../../{input.contigs} ../../../../../{input.reference} && "
        "cd ../../../../../ && mv {params.outdir}/ragoo_output/ragoo.fasta {output}) > {log} 2>&1"


rule filter_chr0_assembly_comparison:
    input:
        "results/{date}/assembly/{sample}/{assembler}/{sample}.ordered.contigs.fasta",
    output:
        "results/{date}/assembly/{sample}/{assembler}/{sample}.contigs.ordered.filtered.fasta",
    log:
        "logs/{date}/ragoo/{assembler}/{sample}_cleaned.log",
    params:
        sample=lambda wildcards: wildcards.sample,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/ragoo-remove-chr0.py"


rule align_contigs_assembly_comparison:
    input:
        target="resources/genomes/main.fasta",
        query="results/{date}/assembly/{sample}/{assembler}/{sample}.contigs.fasta",
    output:
        "results/{date}/assembly/{sample}/{assembler}/main_{sample}.bam",
    log:
        "results/{date}/assembly/{sample}/{assembler}/main_{sample}.log",
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -ax asm5 {input.target} {input.query} -o {output} 2> {log}"


rule quast_assembly_comparison:
    input:
        fasta="results/{date}/assembly/{sample}/{assembler}/{sample}.contigs.fasta",
        bam="results/{date}/assembly/{sample}/{assembler}/main_{sample}.bam",
        reference="resources/genomes/main.fasta",
    output:
        "results/{date}/assembly/{sample}/{assembler}/quast/report.tsv",
    params:
        outdir=lambda x, output: os.path.dirname(output[0]),
    log:
        "logs/{date}/assembly/quast/{assembler}/{sample}.log",
    conda:
        "../envs/quast.yaml"
    threads: 8
    shell:
        "quast.py --min-contig 1 --threads {threads} -o {params.outdir} -r {input.reference} --bam {input.bam} {input.fasta} > {log} 2>&1"


rule plot_assemblies:
    input:
        initial=lambda wildcards: expand(
            "results/{{date}}/assembly/{sample}/{assembler}/{sample}.contigs.fasta",
            sample=get_samples_for_date(wildcards.date),
            assembler=[
                # "abyss",
                "megahit-std",
                "megahit-meta-large",
                "megahit-meta-sensitive",
                "trinity",
                "velvet",
                "metaspades",
                "coronaspades",
                "spades",
                "rnaviralspades",
            ],
        ),
        final=lambda wildcards: expand(
            "results/{{date}}/assembly/{sample}/{assembler}/{sample}.contigs.ordered.filtered.fasta",
            sample=get_samples_for_date(wildcards.date),
            assembler=[
                # "abyss",
                "megahit-std",
                "megahit-meta-large",
                "megahit-meta-sensitive",
                "trinity",
                "velvet",
                "metaspades",
                "coronaspades",
                "spades",
                "rnaviralspades",
            ],
        ),
        quast=lambda wildcards: expand(
            "results/{{date}}/assembly/{sample}/{assembler}/quast/transposed_report.tsv",
            sample=get_samples_for_date(wildcards.date),
            assembler=[
                # "abyss",
                "megahit-std",
                "megahit-meta-large",
                "megahit-meta-sensitive",
                "trinity",
                "velvet",
                "metaspades",
                "coronaspades",
                "spades",
                "rnaviralspades",
            ],
        ),
    output:
        "results/{date}/plots/all_assemblies_largest_contigs.svg",
        "results/{date}/plots/all_assemblies_table.csv",
        "results/{date}/plots/all_assemblies_N50.svg",
        "results/{date}/plots/all_assemblies_genome_fraction.svg",
    log:
        "logs/{date}/all_assemblies_plot.log",
    params:
        samples=lambda wildcards: get_samples_for_date(wildcards.date),
        assembler=[
            # "abyss",
            "megahit-std",
            "megahit-meta-large",
            "megahit-meta-sensitive",
            "trinity",
            "velvet",
            "metaspades",
            "coronaspades",
            "spades",
            "rnaviralspades",
        ],
        amplicon_state=lambda wildcards: get_list_of_amplicon_states(wildcards.date),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot-assembly-comparison.py"
