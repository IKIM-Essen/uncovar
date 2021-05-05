# rule clip_primer:
#     input:
#         bam=expand(
#             "results/{{date}}/mapped/ref~{ref}/{{sample}}.bam",
#             ref=config["adapters"]["amplicon-reference"],
#         ),
#         bed=config["adapters"]["amplicon-primers"],
#         ref_fasta="resources/genomes/{reference}.fasta".format(
#             reference=config["adapters"]["amplicon-reference"]
#         ),
#     output:
#         sortbam=temp("results/{date}/con-clipped-reads/{sample}.bam"),
#         sortindex=temp("results/{date}/con-clipped-reads/{sample}.bam.bai"),
#         clippedbam=temp("results/{date}/con-clipped-reads/{sample}.primerclipped.bam"),
#         hardclippedbam=temp(
#             "results/{date}/con-clipped-reads/{sample}.primerclipped.hard.bam"
#         ),
#         sorthardclippedbam=temp(
#             "results/{date}/con-clipped-reads/{sample}.primerclipped.hard.sorted.bam"
#         ),
#         fq1="results/{date}/con-clipped-reads/{sample}.1.fastq.gz",
#         fq2="results/{date}/con-clipped-reads/{sample}.2.fastq.gz",
#     log:
#         "logs/{date}/primer-clipping/{sample}.log",
#     params:
#         dir=lambda w, output: os.path.dirname(output.sortbam),
#         bam=lambda w, output: output.sortbam.split("/")[-1],
#         dir_depth=lambda w, output: "".join(
#             ["../"] * (len(output.sortbam.split("/")) - 1)
#         ),
#     conda:
#         "../envs/bamclipper.yaml"
#     threads: 10
#     shell:
#         """
#         samtools sort -@ {threads} -o {output.sortbam} {input.bam} > {log} 2>&1
#         samtools index {output.sortbam} >> {log} 2>&1
#         cd {params.dir}
#         bamclipper.sh -b {params.bam} -p {params.dir_depth}{input.bed} -n {threads} >> {params.dir_depth}{log} 2>&1
#         cd {params.dir_depth}
#         fgbio --sam-validation-stringency=LENIENT ClipBam -i {output.clippedbam} -o {output.hardclippedbam} -H true -r {input.ref_fasta} >> {log} 2>&1
#         samtools sort  -@ {threads} -n {output.hardclippedbam} -o {output.sorthardclippedbam}  >> {log} 2>&1
#         samtools fastq -@ {threads} {output.sorthardclippedbam} -1 {output.fq1} -2 {output.fq2}  >> {log} 2>&1
#         """
