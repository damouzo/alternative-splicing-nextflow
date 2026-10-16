# alternative-splicing-nextflow

Nextflow DSL2 pipeline for differential alternative splicing analysis from bulk RNA-seq. It takes the output of nf-core/rnaseq (BAMs and Salmon quantifications) and produces, for every comparison you define, a set of standardised tables, PDF figures and one self-contained HTML report.

## Tools

| Tool | What it measures | Default |
|------|------------------|---------|
| rMATS-turbo | Event-based splicing (SE, A5SS, A3SS, MXE, RI), annotated and de novo | on |
| MAJIQ V3 | Local splicing variations, Bayesian ΔPSI | on |
| IsoformSwitchAnalyzeR | Isoform switches and their consequences (NMD, ORF, coding potential) | on |
| LeafCutter | Differential intron excision | off |
| PEGASAS | Pathway activity versus PSI | off |
| rmats2sashimiplot | Sashimi plots of the top rMATS events | off |

IsoformSwitchAnalyzeR can optionally add PFAM domains and IUPred3 disordered regions (`--run_isar_full_annotation`).

The report also includes volcano plots, PSI PCA, junction QC, cross-tool gene overlap (UpSet), GO/KEGG enrichment, and a differential expression plus splicing comparison when `--dge_dirs` is given.

For every tool, the effect size is **group2 − group1**: positive values mean higher in the second group of each row of `comparisons.csv`. Multiple-testing correction is done per comparison.

## Requirements

- Nextflow
- Docker or Apptainer
- Output of rnaseq (BAM, BAI, Salmon)
- A MAJIQ academic licence if you run MAJIQ (free, https://majiq.biociphers.org)

## Input

Two CSV files.

`samplesheet.csv`: `sample, condition, replicate, bam, bai, salmon_dir`

`comparisons.csv`: `group1, group2`

## Parameters

Required:

- `--input` samplesheet
- `--comparisons` comparisons file
- `--gtf` annotation
- `--strandedness` `unstranded`, `forward` or `reverse`
- `--read_length`

Also needed in some cases:

- `--genome_fasta` if `--use_gffread true`
- `--majiq_license` or the `MAJIQ_LICENSE_FILE` environment variable if MAJIQ is enabled

Commonly changed:

- `--run_majiq`, `--run_isar`, `--run_leafcutter`, `--run_pegasas`, `--run_sashimi` enable or disable tools
- `--organism` `human` or `mouse`
- `--dge_dirs` differential expression results, for the DE + AS comparison
- `--publish_level` how much of the native tool output is published
- `--report_min_reads` minimum junction reads for ranking tables and PCA

The full list, with defaults and per-tool thresholds, is in [docs/usage.md](docs/usage.md). Templates: [assets/params.yaml](assets/params.yaml) and [demo/params.yaml](demo/params.yaml).

## Running

```bash
# Docker
nextflow run main.nf -profile docker -params-file assets/params.yaml

# HPC with Apptainer
nextflow run main.nf -profile apptainer,slurm -params-file demo/params.yaml

# rMATS only, with sashimi plots
nextflow run main.nf -profile docker -params-file params.yaml \
  --run_majiq false --run_isar false \
  --run_sashimi true --sashimi_top_n 10

# LeafCutter and full ISAR annotation
nextflow run main.nf -profile docker -params-file params.yaml \
  --run_leafcutter true \
  --run_isar_full_annotation true --pfam_hmm /data/Pfam-A.hmm
```

MAJIQ:

```bash
export MAJIQ_LICENSE_FILE=/path/to/majiq_license.key
```

## Output

```
results/
├── deliverables/                    # what you send
│   ├── run_info/                    # sample_index, software_versions, qa_report
│   └── contrasts/<comparison_id>/
│       ├── <id>_splicing_report.html
│       ├── tables/                  # <id>.<tool>.{master,significant,summary}.tsv + cross-tool tables
│       └── plots/                   # one PDF per figure, grouped by tool
└── raw/<tool>/<comparison_id>/      # native tool output, for audit only
```

To send the results:

```bash
mv deliverables <name> && zip -r <name>.zip <name>
```

Column definitions are in [docs/output.md](docs/output.md).

## Containers

Images for isar, report, pegasas and leafcutter are built from `containers/` and published to GHCR by `.github/workflows/build-containers.yml`:

```
ghcr.io/damouzo/alternative-splicing-nextflow/<tool>:latest
```

MAJIQ is not built automatically because of its licence. Build it locally:

```bash
docker build -t your-registry/majiq:3.0 containers/majiq/
```

Any image can be overridden with an environment variable (`ISAR_CONTAINER=...`) or a parameter (`--isar_container ...`).


## Documentation

- [docs/usage.md](docs/usage.md) parameters and usage
- [docs/output.md](docs/output.md) output files
- [CHANGELOG.md](CHANGELOG.md)

## License

MIT