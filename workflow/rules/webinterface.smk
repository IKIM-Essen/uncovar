rule post_report:
    input:
        zip_file=expand("results/{mode}-reports/{{date}}.zip", mode=config["mode"]),
    output:
        "results/{date}/.indicators/webinterface/posted-report.json",
    params:
        url=config["webinterface"]["url"] + "report/",
        token=config["webinterface"]["api-token"],
        project_id={"belongs_to": config["webinterface"]["project-id"]},
        additional_data=lambda w: {
            "description": f"Contains samples: {get_samples_for_date(w.date)}"
        },
    log:
        "logs/results/{date}/post-report.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/webinterface/post-file.py"


rule post_sample:
    output:
        "results/{date}/.indicators/webinterface/posted-sample-{sample}.json",
    params:
        url=config["webinterface"]["url"] + "sample/",
        token=config["webinterface"]["api-token"],
        project_id={"project": config["webinterface"]["project-id"]},
    log:
        "logs/results/{date}/post-sample-{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/webinterface/post-sample.py"


rule post_lineage:
    input:
        pango_call="results/{date}/tables/strain-calls/{sample}.polished.strains.pangolin.csv",
        sample_response="results/{date}/.indicators/webinterface/posted-sample-{sample}.json",
    output:
        "results/{date}/.indicators/webinterface/posted-lineage-{sample}.json",
    params:
        url=config["webinterface"]["url"] + "call/lineage/",
        token=config["webinterface"]["api-token"],
    log:
        "logs/results/{date}/post-sample-{sample}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/webinterface/post-lineage.py"


rule post_results:
    input:
        "results/{date}/.indicators/webinterface/posted-report.json",
        lambda wildcards: expand(
            [
                "results/{{date}}/.indicators/webinterface/posted-sample-{sample}.json",
                "results/{{date}}/.indicators/webinterface/posted-lineage-{sample}.json",
            ],
            sample=get_samples_for_date(wildcards.date),
        ),
    output:
        touch("results/{date}/.indicators/webinterface/posted-all-data.json"),
    log:
        "logs/results/{date}/post-sample.log",
