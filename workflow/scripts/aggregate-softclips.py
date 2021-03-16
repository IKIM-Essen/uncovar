sys.stderr = open(snakemake.log[0], "w")
#parameter = snakemake.params.get("parameter", "")


def sum_up_softclips(sm_input, sm_output):
    softclip_dict = {}
    for input in sm_input:
        with open(input, "r") as input_file:
            for line in input_file:
                seq, count = [ele.strip() for ele in line.split(":")]

                if seq not in softclip_dict.keys():
                    softclip_dict[seq] = int(count)
                else:
                    softclip_dict[seq] += int(count)

    softclip_dict = {k: v for k, v in sorted(softclip_dict.items(), key=lambda item: item[1], reverse=True)}
    with open(sm_output, "w") as out_file:
        [print(f"{seq}: {count}", sep="\n", file=out_file) for seq, count in softclip_dict.items()]


if __name__ == "__main__":
    sum_up_softclips(snakemake.input, snakemake.output[0])