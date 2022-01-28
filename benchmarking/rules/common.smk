def get_fastq_pass_path_barcode(wildcards, sample=None):
    if sample is not None:
        return pep.sample_table.loc[sample]["fastq_pass"]
    return pep.sample_table.loc[wildcards.sample]["fastq_pass"]


# def get_fastq_pass_path(wildcards):
#     return os.path.dirname(get_fastq_pass_path_barcode(wildcards))


def get_fast5_pass_path_barcode(wildcards):
    return pep.sample_table.loc[wildcards.sample]["fast5_pass"]


# def get_fast5_pass_path(wildcards):
#     return os.path.dirname(get_fast5_pass_path_barcode(wildcards))


# def get_fastq_input_folder(tech):
#     dir_names = pep.sample_table.loc[
#         pep.sample_table["technology"] == tech, "fq1"
#     ].apply(lambda x: os.path.dirname(x))
#     dir_names = dir_names.unique()
#     assert len(dir_names) == 1, "Can not process files in different dirs."
#     return dir_names


def get_seq_summary(wildcards):
    return pep.sample_table.loc[wildcards.sample]["seq_summary"]


def get_barcode_for_viralrecon_nanopore_sample(wildcards):
    barcode = os.path.basename(os.path.normpath(get_fastq_pass_path_barcode(wildcards)))
    barcode = barcode.replace("barcode0", "")
    barcode = barcode.replace("barcode", "")
    return f"sample,barcode\n{wildcards.sample},{barcode}"


def get_barcode_for_viralrecon_illumina_sample(wildcards):
    fq1, fq2 = get_fastqs(wildcards)
    return f"sample,fastq_1,fastq_2\n{wildcards.sample},{fq1},{fq2}"


def get_fastq_or_fast5(wildcards):
    if wildcards.folder == "fastq_pass":
        return get_fastq_pass_path_barcode(wildcards)
    if wildcards.folder == "fast5_pass":
        return get_fast5_pass_path_barcode(wildcards)


def get_barcode(wildcards):
    return os.path.basename(os.path.normpath(get_fastq_pass_path_barcode(wildcards)))


def get_barcodes(wildcards):
    fastq_paths = [
        get_fastq_pass_path_barcode(None, sample)
        for sample in get_nanopore_samples(wildcards)
    ]
    return [os.path.basename(os.path.normpath(sample)) for sample in fastq_paths]


def get_covpipe_names(wildcards):
    return [sample.replace("_", "__") for sample in get_illumina_samples(wildcards)]


def get_nanopore_samples(wildcards):
    return pep.sample_table.loc[
        pep.sample_table["technology"] == "ont", "sample_name"
    ].values


def get_illumina_samples(wildcards):
    return pep.sample_table.loc[
        pep.sample_table["technology"] == "illumina", "sample_name"
    ].values
