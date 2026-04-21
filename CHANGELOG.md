# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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

## [Unreleased]

### Changed
- Awaiting first post-release updates

---

**Note**: For detailed usage instructions, see [docs/usage.md](docs/usage.md). For output file descriptions, see [docs/output.md](docs/output.md).
