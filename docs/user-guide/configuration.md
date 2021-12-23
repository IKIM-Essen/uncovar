# Configuration

## Execution Mode

Accepted values: `patient`, `environment`. Defaults to `patient`.

Defines the execution mode of UnCoVar.

When the mode is set to `patient`, the sample is assumed come be from a single
host organism and contains only one strain of SARS-CoV-2. The parts of the
workflow for reconstructing the SARS-CoV-2 strain genome are activated.

If the mode is set to `environment`, the sample is assumed to be from the
environment (e.g. wastewater) and to contain different SARS-CoV-2 strains.
The parts of the workflow responsible for creating and analysing individual
genomes (e.g. assembly, lineage calling via Pangolin) are disabled.
