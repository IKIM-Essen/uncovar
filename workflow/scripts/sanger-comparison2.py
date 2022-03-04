import pysam
import pandas as pd


def get_hgvsp(input):
    ann = input.split("|")
    hgvsp = ann[11]
    feature = ann[3]
    if hgvsp:
        # TODO think about regex instead of splitting
        _prefix_enssast_id, alteration = hgvsp.split(":", 1)
        _prefix, alteration = alteration.split(".", 1)
        for triplet, amino in AA_ALPHABET_TRANSLATION.items():
            alteration = alteration.replace(triplet, amino)

        hgvsp = f"{feature}:{alteration}"
        return hgvsp
    else:
        return


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

data = pd.DataFrame(columns=["Sanger", "NGS", "Position"])


start_stop = []
with pysam.AlignmentFile(snakemake.input.sanger_aln, "r") as alignment:
    for aln in alignment.fetch():
        start_stop.append((aln.reference_start, aln.reference_end))

with pysam.VariantFile(snakemake.input.sanger_calls, "rb") as sanger_calls:
    for record in sanger_calls:
        for ann in record.info["ANN"]:
            if get_hgvsp(ann):
                data = pd.concat(
                    [
                        data,
                        pd.DataFrame(
                            {"Sanger": "x", "Position": record.pos},
                            index=[get_hgvsp(ann)],
                        ),
                    ]
                )

with pysam.VariantFile(snakemake.input.ngs_calls, "rb") as ngs_calls:
    for record in ngs_calls:
        for elem in start_stop:
            if record.pos in range(elem[0], elem[1]):
                for ann in record.info["ANN"]:
                    if get_hgvsp(ann):
                        if get_hgvsp(ann) in data.index:
                            data.loc[get_hgvsp(ann)]["NGS"] = "x"
                        else:
                            data = pd.concat(
                                [
                                    data,
                                    pd.DataFrame(
                                        {"NGS": "x", "Position": record.pos},
                                        index=[get_hgvsp(ann)],
                                    ),
                                ]
                            )

data.index.name = "Mutations"
data.to_csv(snakemake.output[0], index=True)
