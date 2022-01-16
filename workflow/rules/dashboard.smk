rule post_report:
    input:
        zip_file="results/reports/{date}.zip",
    output:
        "results/{date}/.indicators/dashboard-post-report.json",
    params:
        url="",
        token="",
        additonal_data=lambda w: {
            "description": f"Contains samples: {get_samples_for_date(w.date)}"
        },
    log:
        "logs/results/{date}/post-report.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/post-file.py"
