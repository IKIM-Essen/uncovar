# Copyright 2022 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.


rule update_sample:
    input:
        "config/pep/samples.csv",
    log:
        "logs/sample_update/preprocessing/sample_csv_update.txt",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/update-sample-sheet.py"


rule barcode_rename:
    input:
        barcodes="barcode-rename.csv",
        #barcode_dir=config["barcode_dir"],
    output:
        #out_dir=directory(config["output_dir"]),
    log:
        "logs/barcode_rename.txt",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/barcode_rename.py"
