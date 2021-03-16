sys.stderr = open(snakemake.log[0], "w")
# parameter = snakemake.params.get("parameter", "")

import pandas as pd
import altair as alt


def mask(string):
    if string == "other":
        return "other"
    elif string == "unmapped":
        return "unmapped"
    elif string == "B.1.1.7":
        return "B.1.1.7"
    elif string == "B.1.351":
        return "B.1.351"
    else:
        return "some strain"


def plot_error(sm_input, sm_output, type="heatmap"):
    results_df = pd.read_csv(sm_input, delimiter="\t")

    results_df["Kallisto output"] = results_df["target_id"].apply(lambda x: mask(x))
    no_of_mixs = int(results_df["mix"].max() + 1)
    
    if type=="heatmap":
        plot = alt.Chart(results_df).mark_rect().encode(
        alt.X('true_fraction:Q', bin=alt.Bin(maxbins=20)),
        alt.Y('est_fraction:Q', bin=alt.Bin(maxbins=20)),
        alt.Color('count(true_fraction):Q', scale=alt.Scale(scheme='greenblue'))
        )

    else: 
        plot = (
            alt.Chart(results_df)
            .mark_point()
            .encode(x="true_fraction:Q", y="est_fraction:Q", color="Kallisto output:N")
        )

    line = (
        alt.Chart(pd.DataFrame({"true_fraction": [0, 1], "est_fraction": [0, 1]}))
        .mark_line()
        .encode(alt.X("true_fraction"), alt.Y("est_fraction"))
    )

    (plot + line).properties(title=f"No. of mixtures {no_of_mixs}").save(sm_output)


if __name__ == "__main__":
    plot_error(snakemake.input[0], snakemake.output[0])
