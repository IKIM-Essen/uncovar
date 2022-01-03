# Copyright 2021 Thomas Battenfeld, Alexander Thomas, Johannes Köster.
# Licensed under the BSD 2-Clause License (https://opensource.org/licenses/BSD-2-Clause)
# This file may not be copied, modified, or distributed
# except according to those terms.

import sys

sys.stderr = open(snakemake.log[0], "w")

import altair as alt
import pandas as pd
import pysam

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


def get_calls():
    variants = []
    for file, date, sample in zip(
        snakemake.input.bcf, snakemake.params.dates, snakemake.params.samples
    ):
        with pysam.VariantFile(file, "rb") as infile:
            for record in infile:
                vaf = record.samples[0]["AF"][0]
                for ann in record.info["ANN"]:
                    ann = ann.split("|")
                    hgvsp = ann[11]
                    enssast_id = ann[6]
                    feature = ann[3]
                    orf = ann[3]
                    if hgvsp:
                        # TODO think about regex instead of splitting
                        enssast_id, alteration = hgvsp.split(":", 1)
                        _prefix, alteration = alteration.split(".", 1)
                        for triplet, amino in AA_ALPHABET_TRANSLATION.items():
                            alteration = alteration.replace(triplet, amino)

                        variants.append(
                            {
                                "feature": feature,
                                "alteration": alteration,
                                "vaf": vaf,
                                "date": date,
                                "sample": sample,
                                "orf": orf,
                            }
                        )
    variants_df = pd.DataFrame(variants)
    variants_df = variants_df[variants_df["orf"] == snakemake.wildcards.ORFNAME]
    return variants_df


def plot_variants_over_time(sm_output, sm_output_table):
    calls = get_calls()

    # write out as table
    calls.to_csv(sm_output_table)

    if len(calls) > 0:
        # get occurrences
        calls["total occurrence"] = calls.groupby("alteration", as_index=False)[
            "alteration"
        ].transform(lambda s: s.count())

        # mask low occurrences
        print(calls["alteration"].value_counts())
        df = pd.DataFrame(calls["alteration"].value_counts())
        df.sort_values(by=["alteration"])
        if len(df.index) > 10:
            # print(calls.loc[calls["alteration"].isin(df.head(10).index)])
            calls.loc[
                ~calls["alteration"].isin(df.head(10).index), "alteration"
            ] = "other occ."
        else:
            calls.loc[calls["total occurrence"] < 0, "alteration"] = "other occ."

    calls.rename(columns={"alteration": "Alteration", "date": "Date"}, inplace=True)

    area_plot = (
        alt.Chart(calls)
        .mark_bar(opacity=0.8)
        .encode(
            x=alt.X("Date:O"),
            y=alt.Y(
                "count()",
                stack="normalize",
                axis=alt.Axis(format="%"),
                title="Fraction in Run",
            ),
            stroke="Alteration",
            color=alt.Color(
                "Alteration",
                scale=alt.Scale(scheme="tableau10"),
                legend=alt.Legend(orient="top"),
            ),
        )
    ).properties(width=800)

    area_plot.save(sm_output)


plot_variants_over_time(snakemake.output[0], snakemake.output[1])
