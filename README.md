# alternative-splicing-nextflow

Standalone Nextflow DSL2 pipeline for differential alternative splicing analysis.

This repository focuses only on alternative splicing. It does not require direct communication with any upstream pipeline. You provide file paths in a samplesheet, and the pipeline consumes those files.

## What It Runs

- rMATS-turbo for event-based splicing analysis
- MAJIQ for local splicing variation and delta PSI
- IsoformSwitchAnalyzeR for isoform switching and consequence analysis
- One consolidated HTML report per comparison

## Input Contract

The pipeline needs two CSV files:

1. samplesheet.csv
2. comparisons.csv


## Required Run Parameters

- --input
- --comparisons
- --gtf
- --strandedness (unstranded, forward, reverse)
- --read_length

If use_gffread=true, you must also provide --genome_fasta.
If run_majiq=true, you must provide --majiq_license.

## Quick Run

Parameter templates included:

- Generic template: [assets/params.yaml](assets/params.yaml)
- Demo template (ad vs old): [demo/params.yaml](demo/params.yaml)

Run with generic template:

nextflow run main.nf -profile docker -params-file assets/params.yaml

Run with demo template:

nextflow run main.nf -profile docker -params-file demo/params.yaml

Or on HPC with Apptainer/Singularity:

nextflow run main.nf -profile apptainer -params-file demo/params.yaml -c custom.config

## Internal Assets

The folder [assets/empty](assets/empty) is internal pipeline scaffolding.

- [assets/empty/NO_RMATS](assets/empty/NO_RMATS)
- [assets/empty/NO_MAJIQ](assets/empty/NO_MAJIQ)
- [assets/empty/NO_ISAR](assets/empty/NO_ISAR)

These sentinel directories are used when a branch is disabled to keep report input wiring stable. They should not be removed.

## MAJIQ License and Runtime

MAJIQ does not need a user/password in params. It needs a valid MAJIQ license key file.

1. Export the license in your shell.

export MAJIQ_LICENSE_FILE=/path/to/licenses/majiq_license.key

2. Run the pipeline.

nextflow run main.nf -profile apptainer -params-file demo/params.yaml

The pipeline now picks `MAJIQ_LICENSE_FILE` automatically as default `majiq_license`.

About installation on HPC:

- You do not need MAJIQ installed system-wide if you run with containers.
- You do need a working container runtime (Apptainer/Singularity or Docker).
- The container images used by the pipeline must be available to the runtime.

Example local image build names used by this pipeline:

- `local/majiq:3.0`
- `local/isar:2.10.0`

If needed, build them from:

- [containers/majiq/Dockerfile](containers/majiq/Dockerfile)
- [containers/isar/Dockerfile](containers/isar/Dockerfile)

## Output Layout

results/
  rmats/<comparison_id>/
  majiq/<comparison_id>/
  isoformswitchr/<comparison_id>/
  report/<comparison_id>_splicing_report.html


## License

MIT
