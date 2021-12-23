# Configuration

## Adapters

There are three ways to transfer adapter sequences to UnCoVar for removing them
from the raw data.

### Config File

The adapter sequences used can be specified in the config file under
`preprocessing` -> `kit adapters`.

For **paired-end data**, the adapters can be detected by per-read overlap
analysis, which seeks the overlap of each pair of reads. The adapter sequences
can specify the adapter sequences for read one by `—adapter_sequence` and for
read two by`—adapter_sequence_r2`. An example is for [Illuminas TruSeq library] (<https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/truseq-rna-v2.html>)
is shown below:

```yaml
"--adapter_sequence = AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
--adapter_sequence_r2 = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT”
```

Adapters for **single-end data** can be specified only using the
`—adapter_sequence` option.

```yaml
"--adapter_sequence = AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
```

### Sample Sheet

The second way to remove adapter sequences is to specify the adapter sequence
directly sample in the sample sheet. The adapters must enter it in a column
called `adapters`. For paired-end and single-end format, see above. Here is
an exemplary samples sheet:

| sample_name | fq1         | fq2         | date       | is_amplicon_data | technology | adapters                                           |
|-------------|-------------|-------------|------------|------------------|------------|----------------------------------------------------|
| example-1   | PATH/TO/fq1 | PATH/TO/fq2 | 1970-01-01 | 1                | illumina   | --adapter_sequence=ACGT --adapter_sequence_r2=TGCA |
| example-2   | PATH/TO/fq  |             | 1970-01-01 | 1                | ion        | --adapter_sequence=ACGT                            |

If an adapter sequence is entered into the sample sheet for one sample, this
adapter sequence is obviously used to trim the sequences of this sample. For
empty entries, UnCoVar uses the adapter sequence from the config file.

### Pre-Defined Adapters

UnCoVar supports two different sequencing kits and their respective adapters,
namely:

1. [Revelo RNA-Seq library preparation kit](https://lifesciences.tecan.com/revelo-rna-seq-library-prep-kit?p=tab--5)
1. [EasySeq RC-PCR SARS CoV-2  Whole Genome Sequencing kit](https://www.nimagen.com/shop/products/rc-cov096/easyseq-sars-cov-2-novel-coronavirus-whole-genome-sequencing-kit)

The `adapters` column in the sample sheet is used to trim the adapters sequences
of these kits. Revelo adapters are trimmed by specifying
`revelo-rna-seq` in the column per sample, while the Nimagen adapters are
removed by specifying `nimagen-easy-seq`.  A short example:

| sample_name | fq1         | fq2         | date       | is_amplicon_data | technology | adapters         |
|-------------|-------------|-------------|------------|------------------|------------|------------------|
| example-1   | PATH/TO/fq1 | PATH/TO/fq2 | 1970-01-01 | 0                | illumina   | revelo-rna-seq   |
| example-2   | PATH/TO/fq  | PATH/TO/fq2 | 1970-01-01 | 1                | illumina   | nimagen-easy-seq |
