# Benchmark

The UnCoVar repository includes benchmarking rules to compare the following workflows:

- [ARTIC fieldbioinformatics](https://github.com/artic-network/fieldbioinformatics)
- [CoVpipe](https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe)
- [HaVoc](https://bitbucket.org/auto_cov_pipeline/havoc/src/master/)
- [ncov2019-artic-nf](https://github.com/connor-lab/ncov2019-artic-nf)
- [nf-core/viralrecon](https://github.com/nf-core/viralrecon)
- [poreCov](https://github.com/replikation/poreCov)
- [Signal](https://github.com/jaleezyy/covid-19-signal)
- [SnakeLines](https://github.com/jbudis/snakelines)
- [UnCoVar](https://github.com/IKIM-Essen/uncovar)
- [V-pipe](https://github.com/cbg-ethz/V-pipe)

The benchmark consisted of

- a precision/recall analysis of the variant call files,
- an execution time overview,
- lineage call comparison based on pangolin, and finally
- a sequence analysis based on quast metrics of the consensus sequences.

To get the above-mentioned analysis run:

```bash
snakemake --cores all --use-conda --snakefile benchmarking/Snakefile  --resources nextflow=1 signal=1 vpipe=1 snakeline=1 uncovar=1
```

The plots can then be found under `results/benchmarking/plots`.

It can happen that the hard disk space is not sufficient. Here it helps to
delete the "work" folder of nextflow manually to free up a few GB.
