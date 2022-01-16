import sys

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import requests

url = snakemake.params.url
token = snakemake.params.token

pango_call = pd.read_csv(snakemake.input.pango_call, squeeze=True)
sample_response = pd.read_json(
    snakemake.input.sample_response, orient="records", typ="series"
)

headers = {"Authorization": f"token {token}"}

data = {}
data.update({"sample": sample_response.id})
data.update(pango_call.to_dict(orient="records")[0])

print(f"Headers: {headers}", file=sys.stderr)
print(f"Data: {data}", file=sys.stderr)


r = requests.post(url, headers=headers, data=data)
r.raise_for_status()

with open(snakemake.output[0], "w") as f:
    print(f"{r.text}", file=f)
