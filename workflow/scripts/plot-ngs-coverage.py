import os

import altair as alt
import pandas as pd
from numpy import cov
from pandas.io.parsers import count_empty_vals

count_dict = {}


for vars in snakemake.input.variants:
    if os.stat(vars).st_size != 0:
        var_df = pd.read_csv(vars, header=None)
        for index, cols in var_df.iterrows():
            if count_dict.get(cols[4]):
                count_dict[cols[4]][1] += 1
            else:
                count_dict[cols[4]] = [0, 1]
            if bool(cols[3]):
                count_dict[cols[4]][0] += 1
        with open(snakemake.output.table, "w") as outfile:
            print(var_df, file=outfile)

sum_dict = {}
count_true = 0
count_all = 0
for key in sorted(count_dict.keys()):
    print(
        key,
        ":",
        count_dict[key][0] + count_true,
        count_dict[key][0],
        count_dict[key][1] + count_all,
        count_dict[key][1],
    )
    sum_dict[key] = [count_dict[key][0] + count_true, count_dict[key][1] + count_all]
    count_true += count_dict[key][0]
    count_all += count_dict[key][1]

count_df = pd.DataFrame.from_dict(count_dict, orient="index", columns=["true", "all"])
sum_df = pd.DataFrame.from_dict(sum_dict, orient="index")
print(sum_df)

# plot = (alt.Chart(count_df.reset_index()).mark_line().transform_window(
#     # Sort the data chronologically
#     sort=[{'field': 'index'}],
#     # Include all previous records before the current record and none after
#     # (This is the default value so you could skip it and it would still work.)
#     frame=[None, 0],
#     # What to add up as you go
#     trues='sum(true)',
#     alls='sum(all)',
#     percent='datum.sum(true) / datum.sum(all)'
# ).encode(
#     x='index:Q',
#     # Plot the calculated field created by the transformation
#     y=alt.Y('alls:Q', axis=alt.Axis(format='%', title='percentage')),
# ).properties(width=600))

# plot.save(snakemake.output[1])

sum_df = pd.DataFrame.from_dict(sum_dict, orient="index", columns=["true", "all"])
print(sum_df)

plot = (
    alt.Chart(sum_df.reset_index())
    .mark_line()
    .transform_calculate(percent="datum.true / datum.all")
    .encode(
        x=alt.X("index:Q", axis=alt.Axis(title="coverage")),
        # Plot the calculated field created by the transformation
        y=alt.Y(
            "percent:Q",
            axis=alt.Axis(
                format="%", title="percentage of sanger variants present in NGS genome"
            ),
        ),
    )
    .properties(width=1000)
)

plot.save(snakemake.output[1])
