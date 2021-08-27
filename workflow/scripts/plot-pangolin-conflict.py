sys.stderr = open(snakemake.log[0], "w")
MIXTURE_PART_INDICATOR = snakemake.params.get("separator", "")
MIXTURE_PERCENTAGE_INDICATOR = snakemake.params.get("percentage", "")

import pandas as pd
import altair as alt


######################################################
# Currently only supporting mixtures with one strain #
######################################################


def plot_pangolin_conflict(sm_input, sm_output):
    # aggregate pangolin outputs
    all_sampes = pd.DataFrame()
    for input in sm_input:
        # get actual lineage
        _prefix, true_lineage_with_percent = input.split(MIXTURE_PART_INDICATOR)
        true_lineage, percent = true_lineage_with_percent.split(
            MIXTURE_PERCENTAGE_INDICATOR
        )
        true_lineage = true_lineage.replace("-", ".")
        percent = percent.replace(".strains.pangolin.csv", "")

        # create df for one call
        pangolin_output = pd.read_csv(input)
        pangolin_output["true_lineage"] = true_lineage
        pangolin_output["true_lineage_percent"] = percent
        all_sampes = all_sampes.append(pangolin_output, ignore_index=True)

    all_sampes["correct_lineage_assignment"] = (
        all_sampes["lineage"] == all_sampes["true_lineage"]
    )

    # get share of correct and incorrect calls
    print(
        all_sampes["correct_lineage_assignment"].value_counts(normalize=True),
        file=sys.stderr,
    )

    wrongly_assigned = all_sampes[all_sampes["correct_lineage_assignment"] == False]

    # plot
    source = wrongly_assigned.copy()
    source.rename(
        columns={
            "lineage": "Called lineage",
            "true_lineage": "Actual lineage",
            "note": "Note",
        },
        inplace=True,
    )

    scatter_plot = (
        alt.Chart(source)
        .mark_circle()
        .encode(
            alt.X("Actual lineage:O"),
            alt.Y(
                "Called lineage:O",
            ),
            size="count()",
            color=alt.Color("count()", scale=alt.Scale(scheme="tableau20")),
        )
    )

    scatter_plot.save(sm_output[0])
    wrongly_assigned.to_csv(sm_output[1])


if __name__ == "__main__":
    plot_pangolin_conflict(snakemake.input, snakemake.output)
