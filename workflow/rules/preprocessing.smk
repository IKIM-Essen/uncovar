rule update_sample:
    input:
        "config/pep/samples.csv",
    log:
        "logs/sample_update/preprocessing/sample_csv_update.txt",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/update-sample-sheet.py"
