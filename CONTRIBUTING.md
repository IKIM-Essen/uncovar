# How to contribute

We use `snakemake` to manage this project.
If you don't have `snakemake`, you should install it. Have a look at the [official installation instruction](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html). If you already have conda or mamba installed, run

```bash
    conda create -n snakemake -c bioconda snakemake snakefmt
```

Once installed, activate your `snakemake` enviroment, e.g. by runing `conda activate snakemake`.

If you have Snakemake already installed, make sure that the uncompromising Snakemake code formatter `snakefmt` is available inside this environment via

```bash
    conda install -c bioconda snakefmt
```

## Before submitting

Before submitting your code, please do the following steps:

1. Add any changes you want
1. Add tests, if possible, for the new changes
1. Edit documentation if you have changed something significant
1. Run `snakefmt workflow` and `black workflow` to format your changes.

## Other help

You can contribute by spreading the word about this project.
It would also be a considerable contribution to write a short article on using this project.
You can also share your best practices with us.
