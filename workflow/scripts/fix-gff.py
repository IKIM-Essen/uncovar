import pandas as pd

def set_id(info_string:str, format:str) -> str:
    if "ID=" in info_string:
        return info_string
    
    infos = dict(x.split("=") for x in info_string.split(";"))

    info_string = f"ID={format}-{infos['Name']};" +info_string
    return info_string

gff = pd.read_csv(snakemake.input[0], sep="\t", header=None)

info_column= gff.columns[-1]
format_column=gff.columns[2]

gff[info_column] = gff.apply(
    lambda x: set_id(x[info_column], x[format_column]), axis = 1
)
    

gff.to_csv(snakemake.output[0], header=False, index=False, sep="\t")