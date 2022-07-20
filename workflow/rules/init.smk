include: "ref.smk"

rule make_directories:
    output:
        data=directory("data"),
        incoming=directory("../incoming"),
        archive=directory("../archive")
    log:
        "log/make_directories_.log"
    shell:
        "for dir in data/ ../archive/ ../incoming/; do if [ ! -d ""$dir"" ];"
        " then mkdir ""$dir""; fi done"