import sys
import json
import re

import pandas as pd
import pysam

sys.stderr = open(snakemake.log[0], "w")


KRAKEN_FILTER_KRITERIA = "D"


def iter_with_samples(inputfiles):
    return zip(snakemake.params.samples, inputfiles)


data = pd.DataFrame()

# add numbers of raw and trimmed reads
for sample, file in iter_with_samples(snakemake.input.reads_unfiltered):
    with open(file) as infile:
        number_reads = json.load(infile)
    data = data.append(
        {
            "# raw reads": number_reads["summary"]["before_filtering"]["total_reads"],
            "# trimmed reads": number_reads["summary"]["after_filtering"][
                "total_reads"
            ],
            "sample": sample,
        },
        ignore_index=True,
    )
data.set_index("sample", inplace=True)

# add numbers of reads used for assembly
for sample, file in iter_with_samples(snakemake.input.reads_used_for_assembly):
    with open(file) as infile:
        data.loc[sample, "# used reads"] = int(infile.read()) * 2


def register_contig_lengths(assemblies, name):
    for sample, file in iter_with_samples(assemblies):
        with pysam.FastxFile(file) as infile:
            data.loc[sample, name] = max(len(contig.sequence) for contig in infile)


# add lenghts of initial contigs
register_contig_lengths(snakemake.input.initial_contigs, "initial contig (bp)")

# add lenghts of polished contigs
register_contig_lengths(snakemake.input.polished_contigs, "final contig (bp)")

# add kraken estimates
species_columns = pd.DataFrame()
for sample, file in iter_with_samples(snakemake.input.kraken):

    kraken_results = pd.read_csv(
        file,
        delimiter="\t",
        names=["%", "covered", "assigned", "code", "taxonomic ID", "name"],
    )
    kraken_results["name"] = kraken_results["name"].str.strip()

    keep_rows = (
        (kraken_results["code"] == KRAKEN_FILTER_KRITERIA)
        | (
            kraken_results["name"]
            == "Severe acute respiratory syndrome-related coronavirus"
        )
        | (kraken_results["name"] == "unclassified")
    )

    kraken_results = kraken_results.loc[keep_rows, ["%", "name"]].set_index("name").T

    eukaryota = "Eukaryota (%)"
    bacteria = "Bacteria (%)"
    viruses = "Viruses (%)"
    sars_cov2 = "thereof SARS (%)"
    unclassified = "Unclassified (%)"
    colnames = {
        "Eukaryota": eukaryota,
        "Bacteria": bacteria,
        "Viruses": viruses,
        "Severe acute respiratory syndrome-related coronavirus": sars_cov2,
        "unclassified": unclassified,
    }
    kraken_results.rename(columns=colnames, inplace=True)
    kraken_results = kraken_results.reindex(
        columns=[eukaryota, bacteria, viruses, sars_cov2, unclassified]
    ).fillna(0)
    kraken_results["sample"] = sample
    species_columns = species_columns.append(kraken_results, ignore_index=True)

data = data.join(species_columns.set_index("sample"))

# add pangolin results
for sample, file in iter_with_samples(snakemake.input.pangolin):
    pangolin_results = pd.read_csv(file)
    assert (
        pangolin_results.shape[0] == 1
    ), "unexpected number of rows (>1) in pangolin results"
    lineage = pangolin_results.loc[0, "lineage"]
    if lineage == "None":
        pangolin_call = "no strain called"
    else:
        match = re.match(
            "((?P<varcount>\d+/\d+) .+ SNPs$)|(seq_len:\d+)$|($)",
            pangolin_results.fillna("").loc[0, "note"].strip(),
        )
        assert (
            match is not None
        ), "unexpected pangolin note, please update above regular expression"
        varcount = match.group("varcount") or ""
        if varcount:
            varcount = f" ({varcount})"
        pangolin_call = f"{lineage}{varcount}"
    data.loc[sample, "pangolin strain (#SNPs)"] = pangolin_call

# add variant calls
AA_ALPHABET_TRANSLATION = {
    "Gly": "G",
    "Ala": "A",
    "Leu": "L",
    "Met": "M",
    "Phe": "F",
    "Trp": "W",
    "Lys": "K",
    "Gln": "Q",
    "Glu": "E",
    "Ser": "S",
    "Pro": "P",
    "Val": "V",
    "Ile": "I",
    "Cys": "C",
    "Tyr": "Y",
    "His": "H",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Thr": "T",
}

for sample, file in iter_with_samples(snakemake.input.bcf):
    variants_of_interest = []
    other_variants = []
    with pysam.VariantFile(file, "rb") as infile:
        for record in infile:
            vaf = record.samples[0]["AF"][0]
            for ann in record.info["ANN"]:
                ann = ann.split("|")
                hgvsp = ann[11]
                enssast_id = ann[6]
                feature = ann[3]
                if hgvsp:
                    # TODO think about regex instead of splitting
                    enssast_id, alteration = hgvsp.split(":", 1)
                    _prefix, alteration = alteration.split(".", 1)
                    for triplet, amino in AA_ALPHABET_TRANSLATION.items():
                        alteration = alteration.replace(triplet, amino)

                    hgvsp = f"{feature}:{alteration}"
                    entry = f"{hgvsp}:{vaf:.3f}"
                    if alteration in snakemake.params.voc.get(feature, {}):
                        variants_of_interest.append(entry)
                    else:
                        other_variants.append(entry)
    data.loc[sample, "variants of interest"] = " ".join(
        sorted(set(variants_of_interest))
    )
    data.loc[sample, "other variants"] = " ".join(sorted(set(other_variants)))

fmt_dict = {
    int : ["# raw reads", "# trimmed reads", "# used reads", "initial contig (bp)", "final contig (bp)"],
    str : ["pangolin strain (#SNPs)", "variants of interest", "other variants"]
}
for dtype, columns in fmt_dict.items(): 
    data[columns] = data[columns].astype(dtype)

data.to_csv(snakemake.output[0])
