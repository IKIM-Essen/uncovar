import sys
from pydoc import describe

sys.stderr = open(snakemake.log[0], "w")

import pandas as pd
import pysam

ncov2019_info_fields = [
    {
        "id": "REF_DP",
        "number": 1,
        "type": "Integer",
        "description": "",
    },
    {
        "id": "REF_RV",
        "number": 1,
        "type": "Integer",
        "description": "",
    },
    {
        "id": "REF_QUAL",
        "number": 1,
        "type": "Integer",
        "description": "",
    },
    {
        "id": "ALT_DP",
        "number": 1,
        "type": "Integer",
        "description": "",
    },
    {
        "id": "ALT_RV",
        "number": 1,
        "type": "Integer",
        "description": "",
    },
    {
        "id": "ALT_QUAL",
        "number": 1,
        "type": "Integer",
        "description": "",
    },
    {
        "id": "ALT_FREQ",
        "number": 1,
        "type": "Float",
        "description": "",
    },
    {
        "id": "TOTAL_DP",
        "number": 1,
        "type": "Integer",
        "description": "",
    },
    {
        "id": "PVAL",
        "number": 1,
        "type": "Float",
        "description": "",
    },
    {
        "id": "GFF_FEATURE",
        "number": 1,
        "type": "Float",
        "description": "",
    },
    {
        "id": "REF_CODON",
        "number": 1,
        "type": "Float",
        "description": "",
    },
    {
        "id": "REF_AA",
        "number": 1,
        "type": "Float",
        "description": "",
    },
    {
        "id": "ALT_CODON",
        "number": 1,
        "type": "Float",
        "description": "",
    },
    {
        "id": "ALT_AA",
        "number": 1,
        "type": "Float",
        "description": "",
    },
]

variant_calls = pd.read_csv(snakemake.input[0], sep="\t")
variant_calls.rename(columns={"REGION": "CONTIG"}, inplace=True)

header = pysam.VariantHeader()

for infotag in ncov2019_info_fields:
    header.info.add(**infotag)

header.filters.add("FAIL", None, None, "Failed")


for contig in variant_calls["CONTIG"].unique():
    header.contigs.add(contig)


with pysam.VariantFile(snakemake.output[0], "w", header=header) as outvcf:
    for entry, infoblock in zip(
        variant_calls.to_dict("records"),
        variant_calls[
            [
                "REF_DP",
                "REF_RV",
                "REF_QUAL",
                "ALT_DP",
                "ALT_RV",
                "ALT_QUAL",
                "ALT_FREQ",
                "TOTAL_DP",
                "PVAL",
                "GFF_FEATURE",
                "REF_CODON",
                "REF_AA",
                "ALT_CODON",
                "ALT_AA",
            ]
        ].to_dict("records"),
    ):
        if entry["ALT"].startswith("-"):
            replace_string = entry["ALT"].replace("-", entry["REF"])
            alleles = (replace_string, entry["REF"])
        else:
            alleles = (entry["REF"], entry["ALT"])

        record = outvcf.new_record()
        record.contig = entry["CONTIG"]
        record.pos = entry["POS"]
        record.alleles = alleles
        record.filter.add("PASS" if entry["PASS"] == True else "FAIL")
        record.info.update(infoblock)
        outvcf.write(record)


# NC_045512.2	22028	.	GAGTTCA	G	228	.	INDEL;IDV=71;IMF=0.876543;DP=81;ADF=2,53;ADR=2,24;AD=4,77;VDB=0.813734;SGB=-0.693147;MQSB=0.654604;MQ0F=0;AC=1;AN=1;DP4=2,2,53,24;MQ=43	GT:PL:DP:SP:ADF:ADR:AD	1:255,0:81:2:2,53:2,24:4,77
# NC_045512.2	22028	.	G	-AGTTCA	.	PASS	REF_DP=577;REF_RV=189;REF_QUAL=43;ALT_DP=670;ALT_RV=0;ALT_QUAL=20;ALT_FREQ=0.995542;TOTAL_DP=673;PVAL=0;GFF_FEATURE=nan;REF_CODON=nan;REF_AA=nan;ALT_CODON=nan;ALT_AA=nan

# NC_045512.2	28247	.	A	-GATTTC	.	PASS	REF_DP=3041;REF_RV=334;REF_QUAL=60;ALT_DP=4272;ALT_RV=0;ALT_QUAL=20;ALT_FREQ=0.883557;TOTAL_DP=4835;PVAL=0;GFF_FEATURE=nan;REF_CODON=nan;REF_AA=nan;ALT_CODON=nan;ALT_AA=nan
# NC_045512.2	28247	.	AGATTTC	A	228	.	INDEL;IDV=110;IMF=0.916667;DP=120;ADF=12,98;ADR=2,8;AD=14,106;VDB=0.390816;SGB=-0.693147;MQSB=0.893421;MQ0F=0;AC=1;AN=1;DP4=12,2,98,8;MQ=43	GT:PL:DP:SP:ADF:ADR:AD	1:255,0:120:5:12,98:2,8:14,106
