import sys

sys.stderr = open(snakemake.log[0], "w")

import requests

input_file = snakemake.input[0]
url=snakemake.params.url
token = snakemake.params.token
project_id = snakemake.params.project_id
additional_data = snakemake.params.additional_data

headers = {"Authorization": f"token {token}"}
data={}

with open(input_file, "rb") as f:
    files = {"zip_file": f}
    data.update(additional_data)
    data.update(project_id)

    print(f"Uploading {input_file} to {url}", file=sys.stderr)
    print(f"Headers: {headers}", file=sys.stderr)
    print(f"Data: {additional_data}", file=sys.stderr)
    print(f"Files: {files}", file=sys.stderr)

    r = requests.post(url, headers=headers, data=data, files=files)
    r.raise_for_status()

    with open(snakemake.output[0], "w") as f:
        print(f"Code: {r.status_code}", file=f)
        print(f"Response: {r.text}", file=f)


