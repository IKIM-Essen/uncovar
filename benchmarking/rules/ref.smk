rule download_artic_primer_schemes:
    output:
        directory("resources/benchmarking/artic/repo"),
    log:
        "logs/download_artic_primer_schemes.log",
    conda:
        "../envs/git.yaml"
    shell:
        "git clone https://github.com/artic-network/artic-ncov2019.git {output} 2> {log}"


rule download_v_pipe:
    output:
        directory("resources/benchmarking/v-pipe/repo"),
    log:
        "logs/download_v_pipe.log",
    conda:
        "../envs/git.yaml"
    shell:
        "git clone --depth 1 --branch sars-cov2 https://github.com/cbg-ethz/V-pipe.git {output} 2> {log}"


rule download_viralrecon_script:
    output:
        "resources/benchmarking/nf-core-viralrecon/fastq_dir_to_samplesheet.py",
    log:
        "logs/download_viralrecon_script.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(wget -L https://raw.githubusercontent.com/nf-core/viralrecon/master/bin/fastq_dir_to_samplesheet.py -O {output} &&"
        " chmod 755 {output})"
        " > {log} 2>&1"


rule download_ViReflow:
    output:
        "resources/benchmarking/ViReflow/ViReflow.py",
    log:
        "logs/download_ViReflow.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(wget 'https://raw.githubusercontent.com/niemasd/ViReflow/master/ViReflow.py' -O {output} &&"
        " chmod 755 {output})"
        " 2> {log}"


rule download_C_VIEW:
    output:
        "resources/benchmarking/C-VIEW/install.sh",
    log:
        "logs/download_C_VIEW.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(wget 'https://raw.githubusercontent.com/ucsd-ccbb/C-VIEW/main/install.sh' -O {output} &&"
        " chmod 755 {output})"
        " 2> {log}"


rule download_SIGNAL:
    output:
        repo=directory("resources/benchmarking/SIGNAL/repo"),
        snakefile="resources/benchmarking/SIGNAL/repo/Snakefile",
        config="resources/benchmarking/SIGNAL/repo/config.yaml",
        script_dir=directory("resources/benchmarking/SIGNAL/repo/scripts"),
    log:
        "logs/download_SIGNAL.log",
    conda:
        "../envs/git.yaml"
    shell:
        "if [ -d '{output.repo}' ]; then rm -Rf {output.repo}; fi &&"
        "git clone --branch v1.4.4 https://github.com/jaleezyy/covid-19-signal {output.repo} 2> {log}"


rule snakelines_download:
    output:
        repo=directory("resources/benchmarking/snakelines/repo"),
        config="resources/benchmarking/snakelines/repo/example/covseq/covseq.yaml",
        snakeline="resources/benchmarking/snakelines/repo/snakelines.snake",
    log:
        "logs/snakelines_download.log",
    conda:
        "../envs/git.yaml"
    shell:
        "if [ -d '{output.repo}' ]; then rm -Rf {output.repo}; fi &&"
        "git clone https://github.com/thomasbtf/snakelines.git {output.repo} 2> {log}"


rule download_HAVoC:
    output:
        script="resources/benchmarking/havoc/HAVoC.sh",
    log:
        "logs/download_HAVoC.log",
    conda:
        "../envs/git.yaml"
    shell:
        "wget 'https://bitbucket.org/auto_cov_pipeline/havoc/raw/e1b9b2be490a4b51085e0ee5c1b21cf7a2a95b1c/HAVoC.sh' -O {output} 2> {log}"
