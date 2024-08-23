import shutil
import pandas as pd
import gzip
import subprocess
import sys
import os


def get_barcode_dirs(source_directory, barcode_numbers):
    try:
        barcode_dirs = []
        for barcode_number in barcode_numbers:
            # Construct the source file path
            source_file = f"barcode{barcode_number}"
            source_path = os.path.join(source_directory, source_file)

            if os.path.exists(source_path):
                barcode_dirs.append(source_path)
            else:
                print(f"Directory '{source_path}' does not exist.")

        return barcode_dirs

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return []


def concatenate_fastq(bc_directory, outfile):
    input_files = os.path.join(bc_directory, "*.fastq.gz")
    print(input_files)
    out_file = f"{os.path.join(outfile, os.path.split(bc_directory)[1])}_all.fastq"
    print(out_file)
    subprocess.Popen(f"zcat {input_files} > {out_file}", shell=True).wait()


def run_sample_prep(source_directory, barcode_numbers, outfile):
    bc_folder = get_barcode_dirs(source_directory, barcode_numbers)
    for item in bc_folder:
        print(item)
        concatenate_fastq(item, outfile)


def rename_files(final_dir):
    # rename files
    renames = pd.read_csv(barcode_csv)
    renames.reset_index(drop=True, inplace=True)

    rename_dict = dict(zip(renames["barcode"], renames["sample_name"]))

    print(rename_dict)

    files = os.listdir(final_dir)
    for file in files:
        num = file.split("_")[0][-2:]
        print(num)
        print(
            final_dir + file + " " + final_dir + str(rename_dict[int(num)]) + ".fastq"
        )
        os.rename(final_dir + file, final_dir + str(rename_dict[int(num)]) + ".fastq")
    print(files)


config = snakemake.config

barcode_csv = str(snakemake.input.barcodes)
source_path = str(config["source_dir"])
out_dir = str(config["output_dir"])

if not os.path.exists(source_path):
    print(f"Source directory '{source_path}' not found.")

# Check if the destination directory exists, create it if not
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# getting barcode numbers
barcode_csv_ = pd.read_csv(barcode_csv, dtype={"barcode": str})
used_barcodes = barcode_csv_["barcode"]

run_sample_prep(source_path, used_barcodes, out_dir)
rename_files(out_dir)
