rule pangolin:
    input:
        "results/assembly/{sample}_assembly.out/final.contigs.fa",
    output:
        "results/pangolin/{sample}_strain.csv",
    log:
        "logs/pangolin/{sample}.log",
    threads: 8
    conda:
        "../envs/pangolin.yaml"
    shell:
        "pangolin {input} --outfile {output}"
