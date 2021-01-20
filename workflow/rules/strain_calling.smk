rule sourmash_compute:
    input:
        "refs/strains.fa.gz",
    output:
        "sourmash/strains.sig",
    log:
        "logs/sourmash/sourmash-compute.log",
    threads: 2
    params:
        k="31",
        scaled="1000",
        extra="",
    wrapper:
        "v0.69.0/bio/sourmash/compute"


# TODO Alexander and Thomas: add sourmash gather rule
