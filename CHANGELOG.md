# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Deliverables publication layer (reestructure_plan.md, contract v1.0.0):
  - `deliverables/metadata/`: `run_manifest.yaml`, `sample_index.tsv`,
    `tools_matrix.tsv` with the per-tool `known_issues` reliability record.
  - `deliverables/contrasts/<id>/data_tables/`: standardised
    master/significant/summary tables per enabled tool (rMATS, MAJIQ, ISAR,
    LeafCutter, PEGASAS) + `cross_tool.master.tsv` gene-level union.
  - rMATS masters flag the `--cstat` numeric floor (`fdr_floor_flag`) and
    mark `is_novel_splice_site`; significant sets use the same rule as the
    HTML report.
  - `deliverables/contrasts/<id>/plots/sashimi/sashimi_index.tsv`: navigable
    PDF index with outdir-relative paths (fixes PDF discoverability).
  - `deliverables/contrasts/<id>/metadata/contrast_manifest.yaml` listing all
    deliverables, thresholds and known issues per comparison.
  - QA gate: `validate_results_contract.py` (structure + content sanity
    checks: LeafCutter status-filter regression, rMATS FDR==0 fraction, pinned
    ISAR q-values, silent cross-tool dropouts); report at
    `deliverables/metadata/results_contract_report.txt`; structural failures
    abort the run.
- New params: `publish_deliverables`, `publish_level` (`core`|`core_raw`|`full`),
  `publish_raw`, `publish_rmats_jcec`, `publish_rmats_individual_counts`,
  `results_contract_version`, `tool_reliability` (known-issues record).
- Report: new "Data exports" section listing the CORE deliverable paths with
  their reliability caveats, plus an early "Known tool reliability issues"
  summary. `RENDER_REPORT` now emits `(comparison_id, html)`.
- `containers/report/Dockerfile` installs `python3` for the deliverables
  exporters.

### Changed
- rMATS RAW publishing: `*.MATS.JCEC.txt` and `individualCounts.*` are no
  longer published by default (`--publish_rmats_jcec` /
  `--publish_rmats_individual_counts` re-enable them).
- `INPUT_CHECK` emits `samples_full` (meta + bam + bai + salmon_dir) for the
  metadata layer.

### Fixed
- Data-integrity fixes from the internal audit (see `internal_audition.md`):
  rMATS, MAJIQ and LeafCutter no longer rely on two independent `groupTuple`
  calls producing the same internal order (Nextflow only guarantees alignment
  *within* one `groupTuple`). Each tool now groups sample ids, BAMs/juncs and
  conditions in a single tuple so column labels always match the real columns.
  - rMATS: `sample_ids` for the report/PEGASAS derive from the same tuple as
    `b1.txt`/`b2.txt`; `RMATS_POST` additionally writes
    `b1_samples.txt`/`b2_samples.txt` (column → sample mapping) inside the
    published results dir.
  - MAJIQ: `MAJIQ_BUILD` pairs each staged BAM with its sample id via a
    channel-derived map instead of index alignment, so the `.sj` files can
    never be generated with the wrong sample name.
  - LeafCutter: `groups.txt` in `LEAFCUTTER_DS` is built from the real column
    header of the perind counts file (mapping each column name to its
    condition), not from the channel list order.
  - PEGASAS: `prepare_pegasas_inputs.py` now validates the channel-provided
    `--g1-ids`/`--g2-ids` positionally against the materialized
    `b1_samples.txt`/`b2_samples.txt`, aborting on any order mismatch (not
    just cardinality mismatches). PEGASAS now requires those files — they are
    written by the updated `RMATS_POST`, and the pipeline fails early with a
    clear message if they are missing (only affects reuse of rMATS outputs
    produced before this change).
  - Input validation: duplicate BAM basenames within a comparison now abort
    at `INPUT_CHECK` (MAJIQ/rMATS/LeafCutter stage BAMs by basename).
  - rMATS: samples that produce no `.rmats` files are now excluded from the
    BAM lists and sample ids too, keeping `b1.txt`/`b2.txt`,
    `b1_samples.txt`/`b2_samples.txt` and the merged `.rmats` consistent.
  - Report: `resolve_sample_labels` now emits a warning instead of silently
    falling back to `Group1_repN` placeholder labels when the id count does
    not match the rMATS column count.

### Changed
- Removed the **QC Metrics** section from the HTML report and the now-dead
  `--nfcore_multiqc_dir` parameter. QC output (MultiQC) is already delivered by
  the upstream nf-core/rnaseq pipeline.
- Sashimi plots: raised `--sashimi_top_n` default from 10 to 30 events per
  event type (up to 150 plots) since reports embed PNGs rather than PDFs.
- rMATS top-events table now shows the top 100 (was 50) by priority score, and
  removed the "FDR < 1e-300 (suelo del test)" annotation from the volcano.
- rMATS de novo (unannotated) splice-site events are now shown in a dedicated
  searchable table of significant events instead of an inline note.
- MAJIQ and LeafCutter summaries are now rendered as tables instead of plain
  text blocks.
- Cross-tool overlap tables and the shared-gene list are now searchable/filterable
  DataTables.
- Translated remaining non-English warning/notice text in the report to English.
- ISAR switch test now selects the engine per contrast via `--isar_test_method`
  (default `auto`): DEXSeq when the smallest condition has ≤5 replicates, else
  satuRn — mirroring the adequacy rule in ISAR's own `isoformSwitchAnalysisPart1`.
  DEXSeq avoids the locfdr empirical-FDR step in satuRn that was collapsing every
  q-value to a single pinned value (~0.99) in high-isoform runs. For this cohort
  (smallest group = 3 in all 10 contrasts) `auto` resolves to DEXSeq everywhere.
  `--isar_test_method dexseq`/`satuRn` force an engine.

### Added
- DE + AS integration: per-tool source column and a volcano per AS tool.

### Changed
- `rmats_novel_ss` now defaults to `true`: rMATS detects unannotated (de novo)
  splice sites. Relevant for spliceosome mutants (DDX41/DHX34); increases
  runtime ~3-5x. De novo events are flagged as "de novo" in the report and
  excluded from sashimi plots (negative coordinates cannot be drawn).
- The HTML report now embeds sashimi plots as PNGs (`data:image/png` at
  `--sashimi_png_dpi`, default 150) instead of base64 PDFs, making the report
  light and self-contained. Vector PDFs are still published under
  `sashimi/<comparison_id>/`. Falls back to embedded PDFs when `pdftoppm` is
  unavailable (older report container).
- Added automatic quality notices in the report: a low-replication warning
  (n < 4 per group) and an orthogonal-validation disclaimer at the end of each
  tool section.
- PSI PCA now selects the top 2000 events by per-event variance instead of the
  first `head(2000)`, better reflecting global PSI variability.
- Cross-tool overlap now normalises gene symbols (trim/case) and maps symbols
  to unique ENSG IDs via the rMATS geneSymbol→GeneID dictionary, making the
  UpSet/pairwise overlap robust to case and alias mismatches.
- `--report_fdr_cutoff` / `--report_dpsi_cutoff` documented as configurable
  report thresholds.
- `run_sashimi` and `run_leafcutter` now default to `false` in `nextflow.config`.
  They must now be enabled explicitly via `--run_sashimi true` /
  `--run_leafcutter true` or the params file.
- Updated `assets/params.yaml` template to document the optional tool toggles and
  their required inputs (`--salmon_merged_tpm` for PEGASAS, `--de_results` for
  the DE + AS integration).

### Added
- `--sashimi_png_dpi` parameter (default `150`) for the resolution of sashimi
  PNGs embedded in the HTML report.
- `docs/usage.md`: "Experimental design and batch structure" best-practices
  subsection and low-replication note under Common Issues.

### Deferred
- Cross-tool overlap significance testing (`phyper`/`SuperExactTest`),
  event-coordinate matching, effect-size tiers, and SpliceAI/Pangolin are
  tracked in `task.todo` for a future iteration.

## [1.0.0] - 2024-04-09

Initial stable release of alternative-splicing-nextflow - a modular Nextflow pipeline for comprehensive differential alternative splicing analysis from nf-core/rnaseq outputs.

### Added

#### Core Functionality
- **Three independent analytical branches**:
  - rMATS-turbo for annotated AS event detection (SE, A5SS, A3SS, MXE, RI)
  - MAJIQ V3 for Local Splicing Variation (LSV) detection with Bayesian deltaPSI quantification
  - IsoformSwitchAnalyzeR for isoform switch identification with functional consequence annotation

#### Input Validation
- Comprehensive samplesheet validation with file existence checks
- **Critical validation features**:
  - Chromosome naming consistency check (BAM headers vs GTF) to prevent silent analysis failures
  - Read length verification against BAM content using pysam (prevents PSI bias in rMATS/MAJIQ)
  - Biological replicate count warnings (minimum 3 recommended)
  - Comparisons file validation

#### Pipeline Modules
- **rMATS subworkflow**:
  - `RMATS_PREP` - Per-sample junction extraction
  - `RMATS_POST` - Statistical testing with Bayesian model
  - Automatic b1.txt/b2.txt generation for comparison groups
  - Support for novel splice site detection (optional)
  
- **MAJIQ subworkflow**:
  - `MAJIQ_BUILD` - Splice graph construction per sample
  - `MAJIQ_DELTAPSI` - Posterior deltaPSI quantification
  - `MAJIQ_VOILA_TSV` - Export to human-readable TSV format
  - Dynamic configuration file generation with per-sample parameters
  
- **IsoformSwitchAnalyzeR subworkflow**:
  - `ISAR_IMPORT` - Salmon quantification import
  - `ISAR_SWITCH_TEST` - Statistical isoform switch testing (DEXSeq/DRIMSeq)
  - `ISAR_EXTRACT_ORF` - Open reading frame extraction
  - `ISAR_SWITCH_CONSEQUENCES` - Functional consequence prediction (domains, signal peptides, NMD, IDR)
  - Optional gffread integration for transcriptome extraction from genome+GTF

#### Reporting
- **Consolidated R Markdown HTML report** with:
  - Interactive volcano plots (plotly)
  - Sortable/searchable result tables (DT)
  - Distribution visualizations (ggplot2)
  - Cross-tool overlap analysis (UpSetR)
  - Auto-generated methods section with citations
  - Complete session info for reproducibility
- MultiQC integration for QC metrics aggregation

#### Workflow Features
- **Conditional tool execution**: Toggle rMATS/MAJIQ/ISAR independently with `--run_*` parameters
- **Modular design**: Each tool runs independently; failure in one doesn't block others
- Per-comparison analysis with automatic grouping by `comparison_id`
- Graceful handling of missing optional inputs (e.g., Salmon directories)

#### Configuration
- nf-core-style configuration with profiles for Docker, Singularity, Conda
- Comprehensive parameter schema with validation
- Resource labels (low, medium, high) for HPC scheduler integration
- Execution profiles: `test`, `docker`, `singularity`, `conda`

#### Containers
- Custom Docker images:
  - `majiq/Dockerfile` - MAJIQ V3 with Python 3.9 and required dependencies
  - `isar/Dockerfile` - IsoformSwitchAnalyzeR with BiocManager packages
- All containers pinned to specific versions for reproducibility

#### Documentation
- **Comprehensive usage guide** (`docs/usage.md`):
  - Complete parameter reference with descriptions
  - Best practices for reproducibility
  - Execution time estimates
  - Resource requirements
  - Troubleshooting guide with 10+ common scenarios
  
- **Detailed output documentation** (`docs/output.md`):
  - Description of all output files
  - Column-by-column explanation of result tables
  - File format reference
  - Example R/Python analysis snippets
  - Data retention recommendations

- **Demo data** (`data_demo/`):
  - Realistic nf-core/rnaseq output structure
  - 4 BAM files (2 conditions × 2 replicates) from public datasets
  - Synthetic Salmon quantifications with realistic TPM values
  - GRCh38 chr22 reference files (GTF + FASTA)
  - Pre-configured samplesheet and comparisons files
  - Automated setup script (`prepare_and_run.sh`)
  - Helper script for samplesheet generation from nf-core outputs (`bin/generate_samplesheet.py`)

#### Helper Scripts
- `bin/validate_samplesheet.py` - Comprehensive input validation
- `bin/prepare_rmats_input.py` - Generate rMATS b1.txt/b2.txt files
- `bin/prepare_majiq_config.py` - Generate MAJIQ .conf files with correct parameters
- `bin/generate_samplesheet.py` - Auto-generate samplesheet from nf-core/rnaseq outputs
- R scripts for ISAR workflow: `isar_import.R`, `isar_switch_test.R`, `isar_extract_orf.R`, `isar_switch_consequences.R`
- `bin/render_report.Rmd` - Comprehensive R Markdown report template

### Parameters

#### Required
- `--input` - Samplesheet CSV path
- `--comparisons` - Comparisons CSV path
- `--gtf` - Gene annotation GTF file
- `--strandedness` - Library strandedness (unstranded/forward/reverse)
- `--read_length` - Post-trimming read length

#### Optional
- `--genome_fasta` - Reference genome FASTA
- `--transcript_fasta` - Transcript FASTA (alternative to gffread)
- `--use_gffread` - Extract transcriptome using gffread (default: true)
- `--genome_build` - Genome build identifier for MAJIQ (default: 'hg38')
- `--run_rmats` / `--run_majiq` / `--run_isar` - Tool toggles (default: all true)
- `--rmats_novelss` - Enable novel splice site detection (default: false)
- `--majiq_license` - MAJIQ license file path
- `--isar_alpha` - ISAR significance threshold (default: 0.05)
- `--max_cpus` / `--max_memory` / `--max_time` - Resource limits

### Tool Versions

- rMATS-turbo: v4.1.2
- MAJIQ: v2.3
- IsoformSwitchAnalyzeR: v2.0.0
- Nextflow: ≥22.10.0
- MultiQC: v1.14

### Technical Specifications

- **Language**: Nextflow DSL2
- **Minimum Nextflow version**: 22.10.0
- **Container technology**: Docker, Singularity/Apptainer
- **Execution environments**: Local, HPC (SLURM, SGE, PBS), Cloud (AWS, GCP, Azure)
- **Resource management**: Automatic retry with increased resources on failure
- **Caching**: Full Nextflow resume support for failed runs

### Repository Structure

Follows nf-core conventions:
```
alternative-splicing-nextflow/
├── main.nf                    # Entry point
├── nextflow.config            # Main configuration
├── workflows/                 # Top-level workflow
├── subworkflows/local/        # Multi-process workflows
├── modules/local/             # Process definitions
├── bin/                       # Executable scripts
├── conf/                      # Configuration profiles
├── containers/                # Custom Dockerfiles
├── docs/                      # Documentation
├── data_demo/                 # Demo data and examples
└── assets/                    # Static files (samplesheets, schemas)
```

### Known Limitations

- Requires minimum 2 biological replicates per condition (3+ recommended)
- MAJIQ requires academic license (free from https://majiq.biociphers.org)
- IsoformSwitchAnalyzeR external tools (SignalP, DeepTMHMM) require separate licenses
- Novel splice site detection in rMATS significantly increases runtime (3-5×)
- Pipeline assumes coordinate-sorted BAM files from STAR (nf-core/rnaseq output)

### Compatibility

- **Upstream**: Designed to consume nf-core/rnaseq v3.14.0+ outputs
- **Genome builds**: Human (GRCh38/hg38), Mouse (GRCm39/mm10), other vertebrates
- **Annotation formats**: Ensembl GTF, Gencode GTF (must match alignment references)

### Future Enhancements

Planned for v1.1.0:
- Multiple comparison support (more than 2 conditions)
- Integration with leafcutter for intron excision analysis
- Sashimi plot generation for top events
- Gene Ontology enrichment of alternatively spliced genes

---

## [1.1.0] - 2026-05-28

### Added

#### Analysis tools
- **Sashimi plots** via rmats2sashimiplot (`--run_sashimi`): per-event arc plots for top rMATS events; configurable top_n, group labels, and scale factors
- **PEGASAS pathway–splicing correlation** (`--run_pegasas`): Python 3 port of PEGASAS KS enrichment + Pearson correlation of PSI vs pathway activity scores
- **LeafCutter intron excision analysis** (`--run_leafcutter`): regtools junction extraction, leafcutter_cluster_regtools clustering, and leafcutter_ds differential splicing; annotated with gene names from GTF
- **ISAR Tier A full functional annotation** (`--run_isar_full_annotation`): PFAM domain analysis via HMMER/hmmscan (`--pfam_hmm`) and IDR prediction via IUPred3; adds `domains_identified` and `IDR_identified` consequence categories

#### Report
- **PSI PCA** section: per-sample PSI matrix from rMATS IncLevel columns, PCA of splicing profiles
- **Splice junction QC** section: per-sample mean inclusion junction count as coverage proxy
- **Event prioritization score** in rMATS table: `−log10(FDR) × |ΔΨ|`
- **Cross-tool UpSet plot** (UpSetR) including rMATS, MAJIQ, ISAR, and LeafCutter gene sets
- **DE + AS dual-hit volcano** (`--de_results`): integrates DESeq2/edgeR results with AS hits
- **PEGASAS KS score plots** and high-correlation event tables
- **LeafCutter section** with cluster significance table and effect size distribution
- **QC Metrics section** using nf-core/rnaseq MultiQC output (`--nfcore_multiqc_dir`)
- **GO/KEGG enrichment** (clusterProfiler; `--organism` human/mouse)

#### Infrastructure
- **GHCR CI/CD** via `.github/workflows/build-containers.yml`: auto-builds isar, report, pegasas, and leafcutter containers on push to main; manual dispatch per target
- **Multi-comparison support**: samplesheet supports N comparisons in a single run
- **Container override params**: `--isar_container`, `--report_container`, `--pegasas_container`, `--leafcutter_container`
- **MAJIQ_SIF env var** support: `MAJIQ_SIF` environment variable sets container path without explicit param

### Changed
- Default `use_gffread` changed from `true` to `false` (user must opt in)
- Pipeline version bumped to 1.1.0 in manifest
- README rewritten to reflect current feature set and GHCR container strategy

### Fixed
- docs/usage.md: removed non-existent params (`--rmats_cutoff`, `--rmats_min_counts`, `--isar_switch_test_method`); corrected `--rmats_novelss` → `--rmats_novel_ss`
- docs/output.md: added sashimi, PEGASAS, LeafCutter output descriptions and updated directory tree

---

**Note**: For detailed usage instructions, see [docs/usage.md](docs/usage.md). For output file descriptions, see [docs/output.md](docs/output.md).
