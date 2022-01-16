import sys

sys.stderr = open(snakemake.log[0], "w")

import requests

# input_file = snakemake.input
# url=snakemake.params.url
# token = snakemake.params.token
# additional_data = snakemake.params.additional_data

input_file = "results/reports/2021-04-12.zip"
url = ""
token = ""
additional_data = {
    "description": "test api post",
}

headers = {"Authorization": f"token {token}"}

with open(input_file, "rb") as f:
    files = {"zip_file": f}

    print(f"Uploading {input_file} to {url}", file=sys.stderr)
    print(f"Headers: {headers}", file=sys.stderr)
    print(f"Data: {additional_data}", file=sys.stderr)
    print(f"Files: {files}", file=sys.stderr)

    r = requests.post(url, headers=headers, data=additional_data, files=files)

    print(f"Code: {r.status_code}", file=sys.stderr)
    print(f"Response: {r.text}", file=sys.stderr)

    r.raise_for_status()
