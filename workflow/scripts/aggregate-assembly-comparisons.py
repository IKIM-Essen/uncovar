import sys

sys.stderr = open(snakemake.log[0], "w")

from collections import defaultdict
from typing import List

import pandas as pd
import pysam


def aggregate_assembly_comparisons(
    bam_files: List[str], samples: List[str], output: str
):
    data = []
    for sample, bam_file_path in zip(samples, bam_files):
        sample_data = defaultdict()
        with pysam.AlignmentFile(bam_file_path) as bam_file:
            for record in bam_file:
                sample_data["Sample"] = sample
                try:
                    sample_data["Edit distance"] = record.get_tag("NM")
                except KeyError:
                    sample_data["Edit distance"] = "tag 'NM' not present"
                sample_data["Cigarstring"] = record.cigarstring

        data.append(sample_data)

    pd.DataFrame(data).to_csv(output, sep="\t", index=False)


aggregate_assembly_comparisons(
    snakemake.input, snakemake.params.samples, snakemake.output[0]
)
