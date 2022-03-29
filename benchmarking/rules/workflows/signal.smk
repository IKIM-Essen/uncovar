# source: https://github.com/jaleezyy/covid-19-signal


rule download_SIGNAL:
    output:
        repo=directory("resources/benchmarking/SIGNAL/repo"),
        snakefile="resources/benchmarking/SIGNAL/repo/Snakefile",
        config="resources/benchmarking/SIGNAL/repo/config.yaml",
        script_dir=directory("resources/benchmarking/SIGNAL/repo/scripts"),
        update="resources/benchmarking/SIGNAL/repo/scripts/assign_lineages.py",
    log:
        "logs/download_SIGNAL.log",
    conda:
        "../envs/git.yaml"
    shell:
        "if [ -d '{output.repo}' ]; then rm -Rf {output.repo}; fi &&"
        "git clone --branch v1.4.4 https://github.com/jaleezyy/covid-19-signal {output.repo} 2> {log}"


rule SIGNAL_prepare_adapter_file:
    output:
        "resources/benchmarking/SIGNAL/resources/NexteraPE-PE.fa",
    log:
        "logs/SIGNAL_prepare_adapter_file.log",
    conda:
        "../../envs/unix.yaml"
    params:
        adapters=lambda w: get_adapters(w)
        .replace("--adapter_sequence ", ">adapter fwd\n")
        .replace(" --adapter_sequence_r2 ", "\n>adapter rev\n"),
    shell:
        "echo '{params.adapters}' > {output} 2> {log}"


rule SIGNAL_download_gbk:
    output:
        "resources/benchmarking/SIGNAL/resources/MN908947.3.gbk",
    log:
        "logs/SIGNAL_download_gbk.log",
    conda:
        "../../envs/unix.yaml"
    params:
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=MN908947.3&rettype=gb&retmode=txt",
    shell:
        "curl -s '{params.url}' > {output} 2> {log}"


rule SIGNAL_download_gff3:
    output:
        "resources/benchmarking/SIGNAL/resources/MN908947.3.gff3",
    log:
        "logs/SIGNAL_download_gff3.log",
    conda:
        "../../envs/unix.yaml"
    params:
        url="https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=MN908947.3",
    shell:
        "curl -s '{params.url}' > {output} 2> {log}"


rule SIGNAL_download_composite_reference:
    output:
        "resources/benchmarking/SIGNAL/resources/GRC38_no_alt_analysis_set.fna.gz",
    log:
        "logs/SIGNAL_composite_reference.log",
    conda:
        "../../envs/unix.yaml"
    params:
        url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    shell:
        "curl -s '{params.url}' > {output} 2> {log}"


rule SIGNAL_unzip_composite_reference:
    input:
        "resources/benchmarking/SIGNAL/resources/GRC38_no_alt_analysis_set.fna.gz",
    output:
        "resources/benchmarking/SIGNAL/resources/GRC38_no_alt_analysis_set.fna",
    log:
        "logs/SIGNAL_unzip_composite_reference.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "gunzip -dk {input} 2> {log}"


rule SIGNAL_composite_reference:
    input:
        human_genome="resources/benchmarking/SIGNAL/resources/GRC38_no_alt_analysis_set.fna",
        virus_reference="resources/genomes/MN908947.fasta",
    output:
        "resources/benchmarking/SIGNAL/resources/composite_human_viral_reference.fna",
    log:
        "logs/SIGNAL_composite_reference.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "cat {input} > {output} 2> {log}"


rule SIGNAL_bwa_index:
    input:
        "resources/benchmarking/SIGNAL/resources/composite_human_viral_reference.fna",
    output:
        idx=multiext(
            "resources/benchmarking/SIGNAL/resources/composite_human_viral_reference",
            ".amb",
            ".ann",
            ".bwt",
            ".pac",
            ".sa",
        ),
    log:
        "logs/SIGNAL_bwa_index.log",
    benchmark:
        "benchmarks/signal_preprocessing/benchmark.txt"
    params:
        algorithm="bwtsw -b 500000000",
    wrapper:
        "v1.1.0/bio/bwa/index"


# we are using the minikraken 8-GB DB, as the ftp option of kraken2 is
# currently borken. See; # https://github.com/DerrickWood/kraken2/issues/518
rule SIGNAL_link_minikraken2:
    input:
        "resources/minikraken-8GB",
    output:
        directory("resources/benchmarking/SIGNAL/resources/Kraken2/db"),
    log:
        "logs/SIGNAL_minikraken2.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "ln -sr {input} {output}"


# TODO: Use the "orignal" SIGNAL kraken2 database, when kraken2 is fixed,
# rule SINGAL_kraken2_tax:
#     output:
#         directory("resources/benchmarking/SIGNAL/resources/Kraken2/db"),
#     log:
#         "logs/SINGAL_kraken2_tax.log",
#     conda:
#         "../../envs/signal.yaml"
#     shell:
#         "kraken2-build --download-taxonomy --db {output} --threads {threads} --use-ftp 2>{log}"


# rule SINGAL_kraken2_lib:
#     input:
#         "resources/benchmarking/SIGNAL/resources/Kraken2/db",
#     output:
#         touch("resources/benchmarking/SIGNAL/resources/Kraken2/.tax"),
#     log:
#         "logs/SINGAL_kraken2_lib.log",
#     conda:
#         "../../envs/signal.yaml"
#     shell:
#         "kraken2-build --download-library viral --db {input} --threads {threads} --use-ftp 2>{log}"


# rule SINGAL_kraken2_build:
#     input:
#         "resources/benchmarking/SIGNAL/resources/Kraken2/db",
#         "resources/benchmarking/SIGNAL/resources/Kraken2/.tax",
#     output:
#         touch("resources/benchmarking/SIGNAL/resources/Kraken2/.build"),
#     log:
#         "logs/SINGAL_kraken2_build.log",
#     threads: 32
#     conda:
#         "../../envs/signal.yaml"
#     shell:
#         " (kraken2-build --build --threads {threads} --db {input[0]} &&"
#         " kraken2-build --clean --threads {threads} --db {input[0]})"
#         "2>{log}"
#         # "kraken2-build --standard --db {output} --threads {threads}  "


rule SIGNAL_configure_config:
    input:
        config="resources/benchmarking/SIGNAL/repo/config.yaml",
        scheme_bed="resources/nCoV-2019.primer.bed",
        composite_reference="resources/benchmarking/SIGNAL/resources/composite_human_viral_reference.fna",
        composite_reference_idx="resources/benchmarking/SIGNAL/resources/composite_human_viral_reference.amb",
        viral_reference_genome="resources/genomes/MN908947.3.fasta",
        gff3="resources/benchmarking/SIGNAL/resources/MN908947.3.gff3",
        gbk="resources/benchmarking/SIGNAL/resources/MN908947.3.gbk",
        kraken2_db="resources/benchmarking/SIGNAL/resources/Kraken2/db",
        amplicon_loc_bed="resources/nCoV-2019.primer.bed",
    output:
        touch("results/benchmarking/SIGNAL/{sample}/.config-updated"),
    log:
        "logs/SIGNAL_prepare_repo/{sample}.log",
    conda:
        "../../envs/unix.yaml"
    params:
        scheme_bed=lambda w, input: os.path.join("../../../../", input.scheme_bed),
        composite_reference=lambda w, input: os.path.join(
            "../../../../", input.composite_reference.replace(".fna", "")
        ),
        viral_reference_genome=lambda w, input: os.path.join(
            "../../../../", input.viral_reference_genome
        ),
        gff3=lambda w, input: os.path.join("../../../../", input.gff3),
        gbk=lambda w, input: os.path.join("../../../../", input.gbk),
        kraken2_db=lambda w, input: os.path.join("../../../../", input.kraken2_db),
        amplicon_loc_bed=lambda w, input: os.path.join(
            "../../../../", input.amplicon_loc_bed
        ),
    shell:
        "(sed -Ei 's#(^scheme_bed: )(.+)$#\\1{params.scheme_bed}#g' {input.config} &&"
        " sed -Ei 's#(^composite_reference: )(.+)$#\\1{params.composite_reference}#g' {input.config} &&"
        " sed -Ei 's#(^viral_reference_genome: )(.+)$#\\1{params.viral_reference_genome}#g' {input.config} &&"
        " sed -Ei 's#(^viral_reference_feature_coords: )(.+)$#\\1{params.gff3}#g' {input.config} &&"
        " sed -Ei 's#(^breseq_reference: )(.+)$#\\1{params.gbk}#g' {input.config} &&"
        " sed -Ei 's#(^kraken2_db: )(.+)$#\\1{params.kraken2_db}#g' {input.config} &&"
        " sed -Ei 's#(^amplicon_loc_bed: )(.+)$#\\1{params.amplicon_loc_bed}#g' {input.config})"


rule SIGNAL_link_config:
    input:
        config="resources/benchmarking/SIGNAL/repo/config.yaml",
        updated="results/benchmarking/SIGNAL/{sample}/.config-updated",
    output:
        "results/benchmarking/SIGNAL/{sample}/config.yaml",
    log:
        "logs/SIGNAL_prepare_repo/{sample}.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "ln -sr {input.config} {output[0]} 2> {log}"


rule SIGNAL_link_script:
    input:
        "resources/benchmarking/SIGNAL/repo/scripts",
    output:
        directory("results/benchmarking/SIGNAL/{sample}/scripts"),
    log:
        "logs/SIGNAL_prepare_repo/{sample}.log",
    conda:
        "../../envs/unix.yaml"
    shell:
        "ln -sr {input} {output} 2> {log}"


rule SIGNAL_create_sample_sheet:
    input:
        get_fastqs,
    output:
        "results/benchmarking/SIGNAL/{sample}/sample_table.csv",
    log:
        "logs/SIGNAL_create_sample_sheet/{sample}.log",
    params:
        abs_path_fq1=lambda w, input: os.path.join(os.getcwd(), input[0]),
        abs_path_fq2=lambda w, input: os.path.join(os.getcwd(), input[1]),
    conda:
        "../../envs/unix.yaml"
    shell:
        "echo 'sample,r1_path,r2_path\n{wildcards.sample},{params.abs_path_fq1},{params.abs_path_fq2}' > {output} 2> {log}"


rule SIGNAL:
    input:
        snakefile="resources/benchmarking/SIGNAL/repo/Snakefile",
        config="results/benchmarking/SIGNAL/{sample}/config.yaml",
        sample_table="results/benchmarking/SIGNAL/{sample}/sample_table.csv",
        script="results/benchmarking/SIGNAL/{sample}/scripts",
    output:
        outdir=temp(directory("results/benchmarking/SIGNAL/{sample}/results_dir")),
        pangolin="results/benchmarking/SIGNAL/{sample}/results_dir/lineage_assignments.tsv",
        consensus="results/benchmarking/SIGNAL/{sample}/results_dir/all_genomes.fa",
        freebayes_pangolin="results/benchmarking/SIGNAL/{sample}/results_dir/freebayes_lineage_assignments.tsv",
        freebayes_consensus="results/benchmarking/SIGNAL/{sample}/results_dir/all_freebayes_genomes.fa",
        vcf="results/benchmarking/SIGNAL/{sample}/results_dir/{sample}/freebayes/{sample}.variants.norm.vcf",
    log:
        "logs/SIGNAL/{sample}.log",
    conda:
        "../../envs/signal.yaml"
    benchmark:
        "benchmarks/signal/{sample}.benchmark.txt"
    threads: 4
    params:
        out_dir=lambda w, input: os.path.dirname(input.config),
        config=lambda w, input: os.path.basename(input.config),
        cwd=os.getcwd(),
    resources:
        signal=1,
    shell:
        "(cd {params.out_dir} && "
        " snakemake -kp --use-conda "
        " --snakefile {params.cwd}/{input.snakefile}"
        " --configfile {params.config}"
        " --rerun-incomplete"
        " --cores={threads}"
        " --conda-prefix=$PWD/.snakemake/conda all)"
        "> {log} 2>&1"
