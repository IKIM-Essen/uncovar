






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
