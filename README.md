# alternative-splicing-nextflow

Nextflow DSL2 pipeline for differential alternative splicing analysis.

Accepts outputs from nf-core/rnaseq (BAMs + Salmon quantifications). Produces one self-contained HTML report per comparison.

## What It Runs

**Core tools** (all enabled by default):

- rMATS-turbo — event-based splicing (SE, A5SS, A3SS, MXE, RI)
- MAJIQ V3 — Local Splicing Variations with Bayesian deltaPSI
- IsoformSwitchAnalyzeR — isoform switches with functional consequences (NMD, ORF, coding potential)

**Optional tools** (disabled by default):

- Sashimi plots (`--run_sashimi`) — rmats2sashimiplot for top rMATS events
- PEGASAS (`--run_pegasas`) — pathway activity × PSI correlation
- LeafCutter (`--run_leafcutter`) — intron excision differential splicing
- ISAR full annotation (`--run_isar_full_annotation`) — PFAM domains + IUPred3 IDR prediction (Tier A)

**Report features** (always included if data available):

- Interactive volcano plots (plotly)
- Cross-tool gene overlap (UpSetR)
- PSI PCA and splice junction QC
- GO/KEGG enrichment (clusterProfiler)
- DE + AS dual-hit volcano (`--de_results`)

## Input

Two CSV files:

1. `samplesheet.csv` — sample, condition, replicate, bam, bai, salmon_dir
2. `comparisons.csv` — group1, group2

## Required Parameters

- `--input`
- `--comparisons`
- `--gtf`
- `--strandedness` (unstranded, forward, reverse)
- `--read_length`

If `use_gffread=true` (default: false), also provide `--genome_fasta`.
If `run_majiq=true`, provide `--majiq_license` or set `MAJIQ_LICENSE_FILE` env var.

## Quick Run

```bash
# Using a params file (recommended)
nextflow run main.nf -profile docker -params-file assets/params.yaml

# HPC with Apptainer
nextflow run main.nf -profile apptainer,slurm -params-file demo/params.yaml

# rMATS + sashimi plots only
nextflow run main.nf -profile docker \
  --run_majiq false --run_isar false \
  --run_sashimi true --sashimi_top_n 10 \
  -params-file params.yaml

# With LeafCutter and ISAR full annotation
nextflow run main.nf -profile docker \
  --run_leafcutter true \
  --run_isar_full_annotation true --pfam_hmm /data/Pfam-A.hmm \
  -params-file params.yaml
```

Parameter templates:

- Generic: [assets/params.yaml](assets/params.yaml)
- Demo (ad vs old): [demo/params.yaml](demo/params.yaml)

## Containers

Container images are built automatically from `containers/` and published to GitHub Container Registry (GHCR) via `.github/workflows/build-containers.yml` on push to `main`.

| Image | Registry tag | Built from |
|-------|-------------|------------|
| isar | `ghcr.io/damouzo/alternative-splicing-nextflow/isar:latest` | `containers/isar/` |
| report | `ghcr.io/damouzo/alternative-splicing-nextflow/report:latest` | `containers/report/` |
| pegasas | `ghcr.io/damouzo/alternative-splicing-nextflow/pegasas:latest` | `containers/pegasas/` |
| leafcutter | `ghcr.io/damouzo/alternative-splicing-nextflow/leafcutter:latest` | `containers/leafcutter/` |

MAJIQ is excluded from auto-build (requires academic licence). Build locally:

```bash
docker build -t your-registry/majiq:3.0 containers/majiq/
```

Override any image via env var or param:

```bash
export ISAR_CONTAINER=my-registry/isar:custom
# or --isar_container my-registry/isar:custom
```

## MAJIQ License

MAJIQ requires an academic licence (free from https://majiq.biociphers.org).

```bash
export MAJIQ_LICENSE_FILE=/path/to/majiq_license.key
nextflow run main.nf -profile apptainer -params-file demo/params.yaml
```

## Output Layout

```
results/
  rmats/<comparison_id>/
  majiq/<comparison_id>/
  isoformswitchr/<comparison_id>/
  leafcutter/<comparison_id>/           # when --run_leafcutter true
  sashimi_plots/<comparison_id>/        # when --run_sashimi true
  pegasas/<comparison_id>/              # when --run_pegasas true
  report/<comparison_id>_splicing_report.html
```

## Internal Assets

`assets/empty/` holds sentinel directories used internally when a tool branch is disabled. They must not be removed.

- `assets/empty/NO_RMATS`
- `assets/empty/NO_MAJIQ`
- `assets/empty/NO_ISAR`
- `assets/empty/NO_SASHIMI`
- `assets/empty/NO_PEGASAS`
- `assets/empty/NO_LEAFCUTTER`

## Documentation

- Parameters and usage: [docs/usage.md](docs/usage.md)
- Output files: [docs/output.md](docs/output.md)

## License

MIT
