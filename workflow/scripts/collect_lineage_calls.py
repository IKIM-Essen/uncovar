sys.stderr = open(snakemake.log[0], "w")

import pandas as pd

def collect_calls(sm_input, sm_output, states, lineage, number, length):
    # agg pangolin calls
    all_pangolin_calls = [] 

    for file, state in zip(sm_input.pangolin, states):
        call = pd.read_csv(file)
        call["actual_lineage"] = lineage
        call["state"] = state
        all_pangolin_calls.append(call)
    
    pangolin_calls = pd.concat(all_pangolin_calls, axis=0, ignore_index=True)
    pangolin_calls["read_number"] = number
    pangolin_calls["read_length"] = length
    pangolin_calls["correct_call"] = pangolin_calls["lineage"] == pangolin_calls["actual_lineage"]
    pangolin_calls = pangolin_calls[["lineage", "actual_lineage", "read_number", "read_length", "correct_call", "state"]]

    # add kallisto calls
    call = pd.read_csv(sm_input.kallisto[0], sep="\t")
    call = call.iloc[[call['fraction'].idxmax()]]
    call["actual_lineage"] = lineage
    call["read_number"] = number
    call["read_length"] = length
    call["correct_call"] = call["target_id"] == call["actual_lineage"]
    call["state"] = "read"
    call.rename(columns={"target_id":"lineage"}, inplace=True)

    call = call[["lineage", "actual_lineage", "read_number", "read_length", "correct_call", "state"]]

    # bring them together
    call = pangolin_calls.append(call)

    call.to_csv(sm_output, sep="\t", index=False)

    

if __name__ == "__main__":
    collect_calls(
        snakemake.input, 
        snakemake.output[0],
        snakemake.params.get("states", ""),
        snakemake.wildcards.lineage.replace("-","."),
        snakemake.wildcards.number,
        snakemake.wildcards.length
    )


    
