# Benchmarking & Tests

## Variant Call Test Cases

If you suspect that variant calls are wrong, you can create test cases and
commit them to the [varlociraptor team](https://github.com/varlociraptor/varlociraptor).

To create test cases, you need to provide the same sample **twice** to UnCoVar:
once sequenced on the Illumina platform and once sequenced on the Oxford
Nanopore platform. Add a new column to the `config/pep/samples.csv` named
`test_case`. Add a unique identifier to the twice sequenced sample. Here is
an example:

| sample_name     | fq1         | fq2         | date       | is_amplicon_data | technology | test_case   |
|-----------------|-------------|-------------|------------|------------------|------------|-------------|
| illumina-sample | PATH/TO/fq1 | PATH/TO/fq2 | 1970-01-01 | 0                | illumina   | test-case-1 |
| nanopore-sample | PATH/TO/fq  |             | 1970-01-01 | 1                | ont        | test-case-1 |

To generate the test cases, request:

```bash
snakemake --cores all --use-conda generate_test_cases
```

After finished execution, you can find the test cases in `results/testcases`
directory together with a summary file.
