import sys

sys.stderr = open(snakemake.log[0], "w")


with open(snakemake.output.accepted_vcfs, "a") as accepted, open(
    snakemake.output.empty_vcfs, "a"
) as denied:
    for happy_path, vc_record_counts_path in zip(
        snakemake.params.happy_paths, snakemake.input.no_of_records
    ):
        with open(vc_record_counts_path, "r") as count_file:
            if int(next(count_file)) > 0:
                print(happy_path, file=accepted)
            else:
                print(vc_record_counts_path, file=denied)
