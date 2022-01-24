# source: https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html
rule artic_guppyplex:
    input:
        get_fastq_pass_path,
    output:
        "results/benchmarking/artic/guppyplex/{sample}.fasta",
    log:
        "logs/artic_guppyplex/{sample}.log",
    conda:
        "../../envs/artic.yaml"
    shell:
        "artic guppyplex --min-length 400 --max-length 700"
        " --directory {input} --output {output}"
        " > {log} 2>&1"


rule artic_minion:
    input:
        fast5=lambda wildcards: get_fast5_pass_path(wildcards),
        fasta="results/benchmarking/artic/guppyplex/{sample}.fasta",
        repo="resources/benchmarking/artic/repo",
    output:
        vcf="results/benchmarking/artic/minion/{sample}/{sample}.merged.vcf",
        consensus="results/benchmarking/artic/minion/{sample}/{sample}.consensus.fasta",
    log:
        "logs/artic_minion/{sample}.log",
    threads: 16
    conda:
        "../../envs/artic.yaml"
    params:
        primer_schemes=lambda w, input: os.path.join(input.repo, "primer_schemes"),
        medaka_model=config["assembly"]["oxford nanopore"]["medaka_model"],
        outdir=get_output_dir,
        cwd=lambda w: os.getcwd(),
    shell:
        "(cd {params.outdir} &&"
        " artic minion --medaka --medaka-model {params.medaka_model}"
        " --normalise 200 --threads {threads}"
        " --scheme-directory {params.cwd}/{params.primer_schemes}"
        " --read-file {params.cwd}/{input.fasta}"
        " nCoV-2019/V3 {wildcards.sample})"
        " > {params.cwd}/{log} 2>&1"


# TODO hand over all samples at once
# source: https://cbg-ethz.github.io/V-pipe/tutorial/sars-cov2/
rule v_pipe_work_dir:
    input:
        "resources/benchmarking/v-pipe/repo",
    output:
        touch("results/benchmarking/v-pipe/{sample}/.work-dir-created"),
    log:
        "logs/v_pipe_work_dir/{sample}.log",
    conda:
        "../../envs/v-pipe.yaml"
    params:
        v_pipe_copy="results/benchmarking/v-pipe/{sample}/",
        work_dir_path="results/benchmarking/v-pipe/{sample}/work",
    shell:
        "(cp -r {input}/* {params.v_pipe_copy} &&"
        " mkdir -p {params.work_dir_path} &&"
        " cd {params.work_dir_path} &&"
        " ../init_project.sh)"
        " > {log} 2>&1"


rule v_pipe_setup_samples:
    input:
        workdir="results/benchmarking/v-pipe/{sample}/.work-dir-created",
        fastqs=get_fastqs,
    output:
        expand(
            "results/benchmarking/v-pipe/{{sample}}/work/samples/{{sample}}/20200102/raw_data/{{sample}}_R{read}.fastq",
            read=[1, 2],
        ),
    log:
        "logs/v_pipe_setup_samples/{sample}.log",
    conda:
        "../../envs/v-pipe.yaml"
    params:
        fq_dir=lambda w, input: os.path.join(
            os.path.dirname(input.workdir),
            "work",
            "samples",
            w.sample,
            "20200102",
            "raw_data",
        ),
    shell:
        "(mkdir -p {params.fq_dir} &&"
        " gzip -d {input.fastqs[0]} -c > {output[0]} &&"
        " gzip -d {input.fastqs[1]} -c > {output[1]})"
        " 2> {log}"


rule v_pipe_dry_run:
    input:
        expand(
            "results/benchmarking/v-pipe/{{sample}}/work/samples/{{sample}}/20200102/raw_data/{{sample}}_R{read}.fastq",
            read=[1, 2],
        ),
    output:
        "results/benchmarking/v-pipe/{sample}/work/samples.tsv",
    log:
        "logs/v_pipe_dry_run/{sample}.log",
    conda:
        "../../envs/v-pipe.yaml"
    params:
        workdir=lambda w: f"results/benchmarking/v-pipe/{w.sample}/work",
    shell:
        "(cd {params.workdir} &&"
        " ./vpipe --dryrun)"
        " > {log} 2>&1"


rule v_pipe_update_sample_sheet:
    input:
        "results/benchmarking/v-pipe/{sample}/work/samples.tsv",
    output:
        touch("results/benchmarking/v-pipe/{sample}/.edited-sample"),
    log:
        "logs/v_pipe_update_sample_sheet/{sample}.log",
    conda:
        "../../envs/v-pipe.yaml"
    shell:
        "sed -i 's/$/\t150/' {input} 2> {log}"


rule v_pipe_run:
    input:
        "results/benchmarking/v-pipe/{sample}/.edited-sample",
    output:
        vcf="results/benchmarking/v-pipe/{sample}/work/samples/{sample}/20200102/variants/SNVs/snvs.vcf",
        consensus="results/benchmarking/v-pipe/{sample}/work/samples/{sample}/20200102/references/ref_majority.fasta",
    log:
        "logs/v_pipe_run/{sample}.log",
    conda:
        "../../envs/v-pipe.yaml"
    params:
        workdir=lambda w, input: os.path.join(os.path.dirname(input[0]), "work"),
    threads: 8
    shell:
        "(cd {params.workdir} &&"
        " ./vpipe --cores {threads} -p -F)"
        " > {log} 2>&1"


# source: https://github.com/replikation/poreCov
rule poreCov_sample_sheet:
    input:
        get_fastq_pass_path,
    output:
        "results/benchmarking/poreCov/{sample}/sample_names.csv",
    log:
        "logs/poreCov_sample_sheet/{sample}.log",
    conda:
        "../../envs/unix.yaml"
    params:
        barcode=lambda w, input: os.path.basename(input[0]),
    shell:
        "(echo _id,Status,Description > {output} &&"
        " echo {wildcards.sample},{params.barcode},- >> {output})"
        " 2> {log}"


rule poreCov:
    input:
        fastq_pass=get_fastq_pass_path,
        sample_names="results/benchmarking/poreCov/{sample}/sample_names.csv",
        # any --<argname> pipeline file arguments can be given here as <argname>=<path>
    output:
        consensus="results/benchmarking/poreCov/{sample}/2.Genomes/all_consensus_sequences/{sample}.consensus.fasta",
        lineage_call="results/benchmarking/poreCov/28998_ont/3.Lineages_Clades_Mutations/{sample}/lineage_report_{sample}.csv",
    log:
        "logs/poreCov/{sample}.log",
    threads: 8
    params:
        pipeline="replikation/poreCov",
        revision="1.0.0",
        profile=["local", "docker"],
        flags="--update",
        cores=lambda w, threads: threads,
        output=lambda w: f"results/benchmarking/poreCov/{w.sample}",
        samples=lambda W, input: input.sample_names,
        # any --<argname> pipeline arguments can be given here as <argname>=<value>
    handover: True
    conda:
        "../../envs/nextflow.yaml"
    script:
        "../../scripts/benchmarking/nextflow.py"


# source: https://github.com/connor-lab/ncov2019-artic-nf
# rule ncov2019_artic_nf_illumina_data_prep:
#     input:
#         get_fastqs,
#     output:
#         d=directory("data/ref-data/{sample}"),
#         fq1=temp("data/ref-data/{sample}/{sample}_R1.fastq"),
#         fq2=temp("data/ref-data/{sample}/{sample}_R2.fastq"),
#     log:
#         "logs/ncov2019_artic_nf_illumina_data_prep/{sample}.log",
#     conda:
#         "../../envs/unix.yaml"
#     shell:
#         "(mkdir -p {output.d} &&"
#         " gzip -d {input[0]} -c > {output.fq1} &&"
#         " gzip -d {input[1]} -c > {output.fq2})"
#         " 2> {log}"


# TODO need seq summary
# rule ncov2019_artic_nf_illumina:
#     input:
#         directory="data/ref-data/{sample}",
#         # any --<argname> pipeline file arguments can be given here as <argname>=<path>
#     output:
#         consensus="results/benchmarking/ncov2019_artic_nf/{sample}/{sample}.qc.csv",
#     log:
#         "logs/ncov2019_artic_nf/{sample}.log"
#     threads: 8
#     params:
#         pipeline="connor-lab/ncov2019-artic-nf",
#         profile=["conda"],
#         flags="--illumina",
#         outdir = lambda w: f"results/benchmarking/ncov2019_artic_nf/{w.sample}",
#         prefix = lambda w: w.sample,
#         # any --<argname> pipeline arguments can be given here as <argname>=<value>
#     handover: True
#     conda:
#         "../../envs/nextflow.yaml"
#     script:
#         "../../scripts/benchmarking/nextflow.py"

# TODO need seq summary
# use rule ncov2019_artic_nf_illumina as bncov2019_artic_nf_ont with:
#     input:
#         basecalled_fastq=""
#         fast5_pass=""
#     output:
#         consensus="results/benchmarking/ncov2019_artic_nf/nanopore/{sample}/{sample}.qc.csv",
#     log:
#         "logs/ncov2019_artic_nf/nanopore/{sample}.log"
#     params:
#         pipeline="connor-lab/ncov2019-artic-nf",
#         profile=["conda"],
#         flags="--illumina",
#         outdir = lambda w: f"results/benchmarking/ncov2019_artic_nf/illumina/{w.sample}",
#         prefix = lambda w: w.sample,

# source: https://github.com/nf-core/viralrecon
# rule nf_core_viralrecon_sample_sheet:
#     input:
#         script="resources/benchmarking/nf-core-viralrecon/fastq_dir_to_samplesheet.py",
#         fastq_dir="data/fastq_pass",
#     output:
#         "results/benchmarking/nf-core-viralrecon/nanopore/samplesheet.csv",
#     log:
#         "logs/nf_core_viralrecon_sample_sheet.log"
#     conda:
#         "../../envs/python.yaml"
#     shell:
#         "./{input.script} {input.fastq_dir} {output}"

# TODO need seq summary
# rule nf_core_viralrecon:
#     input:
#         input="results/benchmarking/nf-core-viralrecon/samplesheet.csv",
#         # any --<argname> pipeline file arguments can be given here as <argname>=<path>
#     output:
#         consensus="results/benchmarking/nf-core-viralrecon/{sample}/{sample}.qc.csv",
#     log:
#         "logs/nf-core-viralrecon/{sample}.log"
#     threads: 8
#     params:
#         pipeline="nf-core/viralrecon",
#         profile=["docker"],
#         platform="illumina",
#         protocol="metagenomic",
#         genome="'MN908947.3'",
#         # any --<argname> pipeline arguments can be given here as <argname>=<value>
#     handover: True
#     conda:
#         "../../envs/nextflow.yaml"
#     script:
#         "../../scripts/benchmarking/nextflow.py"

# TODO need seq summary
# use rule nf_core_viralrecon as nf_core_viralrecon_nanopore with:
#     input:
#         input="sampleeet.csv",
#         sequencing_summary="",
#         fastq_dir="",
#         fast5_dir="",
#         # any --<argname> pipeline file arguments can be given here as <argname>=<path>
#     output:
#         consensus="results/benchmarking/nf-core-viralrecon/nanopore",
#     log:
#         "logs/nf-core-nf_core_viralrecon_nanopore/{sample}.log"
#     params:
#         pipeline="nf-core/viralrecon",
#         profile=["docker"],
#         platform="nanopore",
#         genome="'MN908947.3'",
#         # any --<argname> pipeline arguments can be given here as <argname>=<value>


# source: https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe#3-usage
rule CovPipe:
    input:
        input_dir=get_fastq_input_folder(ILLUMINA),
        reference="resources/genomes/main.fasta",
    output:
        directory("results/benchmarking/CovPipe"),
    log:
        "logs/CovPipe.log",
    conda:
        "../../envs/covpipe.yaml"
    shell:
        "ncov_minipipe --reference {input.reference}"
        "--input {input.input_dir} "
        "-o {output}"


# TODO Need s3 bucket
# source: https://github.com/niemasd/ViReflow
# rule ViReflow:
#     input:
#         fq=get_fastqs,
#         script="resources/benchmarking/ViReflow/ViReflow.py",
#         reference="resources/genomes/main.fasta",
#         gff="resources/annotation.gff.gz",
#         bed="resources/primers.bed",
#     output:
#         directory("results/benchmarking/ViReflow/{sample}"),
#     log:
#         "logs/ViReflow/{sample}.log",
#     threads: 8
#     conda:
#         "../../envs/python.yaml"
#     shell:
#         "./{input.script} --destination {output} --reference_fasta {input.reference} "
#         "--reference_gff {input.gff} --primer_bed {input.bed} --output {output} "
#         "--threads {threads} --optional_pangolin true  {input.fq[0]}"
# C-View is not installable or runable.
# -> Needs sudo
# -> paths to softwaredir and anaconda dir sometimes hardcoded
# -> Was not able to start
# rule C_VIEW_install:
#     input:
#         "resources/benchmarking/C-VIEW/install.sh",
#     output:
#         directory("resources/benchmarking/C-VIEW/softwaredir"),
#     log:
#         "logs/C_VIEW_install.log"
#     conda:
#         "../../envs/c-view.yaml"
#     shell:
#         "./{input} {output} $(which anaconda | sed 's/\/bin.*//g')"
# What to compare when benchmarking UnCoVar to other pipelines?
# -> de novo sequences, consensus sequences, variants calls, lineage (pangolin) calls of uncovar vs other pipeline.
# How to compare de novo / consensus sequences?
# -> Alignments
# --> Edit distances to SARS-CoV-2 reference genome by uncovar vs. other pipeline
# --> Identity to SARS-CoV-2 reference genome by uncovar vs. other pipeline
# --> Share N in de novo / consensus sequence by uncovar vs. other pipeline
# --> Visualise difference in seqs. by uncovar vs. other pipeline (b.c. of masking in uncovar)
# How to compare variant calls?
# -> Number of variants, which are also present in sanger sequencing
# -> Number of total found variants by uncovar vs. other pipeline
# -> Number of variants that were found by both pipelines
# -> Uniquely found variants by uncovar vs. other pipeline
# -> Are there variants that are unique to certain vocs, which are only found by one pipeline?
# -> How to integrate probs/ vafs / depth other VCF metrics?
# -> Precision und Recall auf den von Sanger abgedeckten bereichen berechnen
# How to compare Pangolin Calls
# -> Pangolin call by uncovar vs. other pipeline
# -> How to compare Kallisto?
