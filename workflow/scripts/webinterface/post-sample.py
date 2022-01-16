import sys

sys.stderr = open(snakemake.log[0], "w")

import requests

url=snakemake.params.url
token = snakemake.params.token
project_id = snakemake.params.project_id
sample=snakemake.wildcards.sample
date=snakemake.wildcards.date

headers = {"Authorization": f"token {token}"}

data={}
data.update(project_id)
data.update({"name" : sample})
data.update({"date" : date})

print(f"Headers: {headers}", file=sys.stderr)
print(f"Data: {data}", file=sys.stderr)

r = requests.post(url, headers=headers, data=data)
r.raise_for_status()

with open(snakemake.output[0], "w") as f:
    print(f"{r.text}", file=f)
