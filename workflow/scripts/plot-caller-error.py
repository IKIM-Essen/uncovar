import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import altair as alt


def mask(strain):
    return (
        strain
        if strain in ["other", "unmapped", "B.1.1.7", "B.1.351"]
        else "some strain"
    )


def plot_error_heatmap(sm_input, sm_output, type="heatmap"):
    results_df = pd.read_csv(sm_input, delimiter="\t")

    results_df["Output"] = results_df["target_id"].apply(lambda x: mask(x))
    results_df = results_df[results_df["Output"] != "unmapped"]
    results_df = results_df[results_df["Output"] != "other"]
    no_of_mixs = int(results_df["mix"].max() + 1)

    if type == "heatmap":
        plot = (
            alt.Chart(results_df)
            .mark_rect()
            .encode(
                alt.X(
                    "true_fraction:Q",
                    axis=alt.Axis(format="%", title="True Fraction"),
                    bin=alt.Bin(maxbins=35),
                    scale=alt.Scale(
                        domain=(0.0, 1.0),
                        bins=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                    ),
                ),
                alt.Y(
                    "est_fraction:Q",
                    axis=alt.Axis(format="%", title="Est. Fraction"),
                    bin=alt.Bin(maxbins=35),
                    scale=alt.Scale(
                        domain=(0.0, 1.0),
                        bins=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
                    ),
                ),
                alt.Color(
                    "count(true_fraction):Q",
                    scale=alt.Scale(type="log", scheme="greenblue"),
                ),
            )
        )

        # Define the degree of the polynomial fits
        degree_list = [1, 3]

        base = (
            alt.Chart(results_df)
            .mark_circle()
            .encode(alt.X("true_fraction"), alt.Y("est_fraction"))
        )

        polynomial_fit = [
            base.transform_regression(
                "true_fraction",
                "est_fraction",
                method="poly",
                order=order,
                as_=["true_fraction", str(order)],
            )
            .mark_line(color="black")
            .transform_fold([str(order)], as_=["Poly. Degree", "est_fraction"])
            .encode(alt.Color("Poly. Degree:N"))
            for order in degree_list
        ]

        # plot = alt.layer(plot, *polynomial_fit)

    else:
        plot = (
            alt.Chart(results_df)
            .mark_point()
            .encode(x="true_fraction:Q", y="est_fraction:Q", color="Kallisto output:N")
        )

    line = (
        alt.Chart(pd.DataFrame({"true_fraction": [0, 1], "est_fraction": [0, 1]}))
        .mark_line(color="grey")
        .encode(alt.X("true_fraction"), alt.Y("est_fraction"))
    )

    (plot + line).properties(
        title=f"{snakemake.wildcards.caller}, No. of mixtures {no_of_mixs}"
    ).save(sm_output)


def plot_bar_of_zeros(sm_input, sm_output):
    results_df = pd.read_csv(sm_input, delimiter="\t")

    results_df["Output"] = results_df["target_id"].apply(lambda x: mask(x))
    results_df = results_df[results_df["Output"] != "unmapped"]
    results_df = results_df[results_df["Output"] != "other"]

    results_df = results_df[
        (results_df["true_fraction"] == 0) | (results_df["est_fraction"] == 0)
    ]

    results_df["Axis"] = results_df["true_fraction"].apply(
        lambda x: "True Fraction = 0" if x == 0 else "Est. Fraction = 0"
    )

    base = alt.Chart(results_df[["true_fraction", "target_id", "Axis"]])

    bar = (
        base.mark_bar()
        .encode(x=alt.X("target_id", sort="-y"), y=alt.Y("count(target_id)"))
        .facet(
            row="Axis:N",
        )
    ).properties(
        title=f"Count of records where the other fraction is 0. Shows the x and y axis of the heatmap in more detail."
    )

    (bar).save(sm_output)


def plot_worst_predictons_content(sm_input, sm_output):
    plots = []

    for i in range(3):
        try:
            results_df = pd.read_csv(sm_input, delimiter="\t")

            results_df["Output"] = results_df["target_id"].apply(lambda x: mask(x))
            results_df = results_df[results_df["Output"] != "unmapped"]
            results_df = results_df[results_df["Output"] != "other"]

            results_df_false = results_df[results_df["true_fraction"] == 0]
            worst_predictions = results_df_false.target_id.value_counts()
            worst_prediction = worst_predictions.index[i]
            worst_predictions = results_df_false[
                results_df_false["target_id"] == worst_prediction
            ].mix.unique()

            results_df = results_df[results_df.mix.isin(worst_predictions)]
            results_df = results_df[results_df.target_id != worst_prediction]
            results_df = results_df[results_df.true_fraction != 0]

            no_samples = len(results_df.mix.unique())

            plot = (
                alt.Chart(results_df)
                .mark_bar()
                .encode(
                    x=alt.X("target_id:O", sort="-y"),
                    y="count(target_id)",
                    color="true_fraction:Q",
                )
                .properties(
                    title=f"Composition of {no_samples} mixtures, where {worst_prediction} is called, but not contained in mixture"
                )
            )

            plots.append(plot)
        except:
            # TODO make exception handling more fine-grained
            pass

    alt.vconcat(*plots).save(sm_output)


if __name__ == "__main__":
    plot_error_heatmap(snakemake.input[0], snakemake.output[0])
    plot_bar_of_zeros(snakemake.input[0], snakemake.output[1])
    plot_worst_predictons_content(snakemake.input[0], snakemake.output[2])
