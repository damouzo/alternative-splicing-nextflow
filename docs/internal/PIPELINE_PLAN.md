# Alternative Splicing Nextflow Pipeline — Design and Implementation Plan

---

## 0. Big-Picture Vision: Modular RNA-seq Analysis Ecosystem

This pipeline is one component of a larger modular analysis strategy. Each major
analytical objective lives in its own independent, versioned Nextflow pipeline that
consumes outputs from the step upstream. This makes each component independently
publishable, re-runnable, and maintainable.

```
FASTQs (Novogene)
        │
        ▼
┌──────────────────────────────────────────┐
│   nf-core/rnaseq  -r 3.14.0              │  ← community-maintained, citable by name+version
│   FastQC → TrimGalore → STAR → Salmon    │
└───────────────┬──────────────────────────┘
                │
       ┌────────┴──────────┬───────────────────────┐
       │                   │                       │
  BAMs + BAI        Salmon quant.sf           (future sources)
       │                   │                       │
       ▼                   ▼                       ▼
┌───────────────┐  ┌────────────────────┐  ┌────────────────────┐
│ alternative-  │  │ deg-gsea-nextflow  │  │pseudotime-nextflow │
│ splicing-     │  │ (future pipeline)  │  │(future pipeline)   │
│ nextflow      │  │ DESeq2 + edgeR     │  │Monocle3 / scVelo   │
│ ← THIS REPO   │  │ fgsea + clusterP.  │  │Palantir            │
└───────────────┘  └────────────────────┘  └────────────────────┘
```

**Why this design?**
- **nf-core/rnaseq** is maintained by the community, benchmarked, and citable in papers
  as "we used nf-core/rnaseq v3.14.0" — reviewers accept this without further
  justification. You never have to maintain the alignment stack.
- **BAMs are produced once** by nf-core/rnaseq and consumed simultaneously by all
  downstream pipelines. No re-alignment, no duplicated compute.
- **Each downstream pipeline has a single responsibility**, making it smaller, easier to
  test, and faster to re-run when only one analysis needs updating.
- **This pipeline never touches raw reads.** It is a pure analysis pipeline, which
  dramatically reduces its complexity and dependency surface.

---

## Overview

This document specifies the architecture, logic, parameter decisions, and implementation
steps for `alternative-splicing-nextflow` — the differential alternative splicing branch
of the modular ecosystem described in Section 0.

**Input**: coordinate-sorted STAR BAMs (+ BAI indices) and Salmon transcript
quantification directories produced by a prior run of nf-core/rnaseq. No alignment is
performed by this pipeline. The pre-existing files are fed directly into three
independent splicing analysis tools.

The three analytical branches are:

1. **rMATS-turbo** — differential alternative splicing across five annotated AS event
   types (SE, A5SS, A3SS, MXE, RI), with optional novel splice-site detection.
2. **MAJIQ V3** — Local Splicing Variation (LSV) detection and deltaPSI quantification
   without assuming identical replicate structures.
3. **IsoformSwitchAnalyzeR** — isoform-switch identification with functional consequence
   annotation (protein domains, signal peptides, intrinsically disordered regions, etc.).

A final R Markdown HTML report consolidates all results.

---

## 1. Input Specification

### 1.1 Input

The pipeline accepts pre-processed outputs from **nf-core/rnaseq** (or any STAR +
Salmon run that produces equivalent files). Required per sample:

| File | Source in nf-core/rnaseq output | Required by |
|------|---------------------------------|-------------|
| Coordinate-sorted BAM | `star_salmon/<sample>/` | rMATS, MAJIQ |
| BAM index (`.bai`) | `star_salmon/<sample>/` | rMATS, MAJIQ |
| Salmon quant directory | `star_salmon/<sample>/` | IsoformSwitchAnalyzeR |

> **nf-core/rnaseq prerequisite**: run with `--aligner star_salmon` (default) and
> `--save_align_intermeds`. Without `--save_align_intermeds`, nf-core/rnaseq deletes
> the STAR BAM files at the end of the run.

### 1.2 Samplesheet Format

A CSV samplesheet is required. Path: passed via `--input`.

```csv
sample,condition,replicate,bam,bai,salmon_dir
ctrl_rep1,control,1,/path/to/rnaseq_results/star_salmon/ctrl_rep1/ctrl_rep1.markdup.sorted.bam,/path/to/rnaseq_results/star_salmon/ctrl_rep1/ctrl_rep1.markdup.sorted.bam.bai,/path/to/rnaseq_results/star_salmon/ctrl_rep1/
ctrl_rep2,control,2,/path/to/rnaseq_results/star_salmon/ctrl_rep2/ctrl_rep2.markdup.sorted.bam,/path/to/rnaseq_results/star_salmon/ctrl_rep2/ctrl_rep2.markdup.sorted.bam.bai,/path/to/rnaseq_results/star_salmon/ctrl_rep2/
ctrl_rep3,control,3,/path/to/rnaseq_results/star_salmon/ctrl_rep3/ctrl_rep3.markdup.sorted.bam,/path/to/rnaseq_results/star_salmon/ctrl_rep3/ctrl_rep3.markdup.sorted.bam.bai,/path/to/rnaseq_results/star_salmon/ctrl_rep3/
treat_rep1,treatment,1,/path/to/rnaseq_results/star_salmon/treat_rep1/treat_rep1.markdup.sorted.bam,/path/to/rnaseq_results/star_salmon/treat_rep1/treat_rep1.markdup.sorted.bam.bai,/path/to/rnaseq_results/star_salmon/treat_rep1/
treat_rep2,treatment,2,/path/to/rnaseq_results/star_salmon/treat_rep2/treat_rep2.markdup.sorted.bam,/path/to/rnaseq_results/star_salmon/treat_rep2/treat_rep2.markdup.sorted.bam.bai,/path/to/rnaseq_results/star_salmon/treat_rep2/
treat_rep3,treatment,3,/path/to/rnaseq_results/star_salmon/treat_rep3/treat_rep3.markdup.sorted.bam,/path/to/rnaseq_results/star_salmon/treat_rep3/treat_rep3.markdup.sorted.bam.bai,/path/to/rnaseq_results/star_salmon/treat_rep3/
```

Rules enforced by the validation script at startup:
- Column names are case-sensitive and must match exactly.
- `sample` values must be unique.
- `bam`, `bai`, and `salmon_dir` must point to existing paths.
- `condition` values must produce exactly two distinct groups (one pairwise comparison).
- Each condition must have a minimum of **3 replicates**. The statistical models of all
  three tools (rMATS Bayesian model, MAJIQ posterior estimation, ISAR satuRn test) require
  at least 3 replicates to produce reliable results. The pipeline will emit a warning if
  either group has fewer than 3 replicates and will exit with an error if either group has
  fewer than 2.

### 1.3 Comparisons File

A two-column CSV passed via `--comparisons`. Specifies which condition is treated as
group 1 (reference/control) and which as group 2 (treatment/condition of interest).

```csv
group1,group2
control,treatment
```

For future multi-group designs this file can hold multiple rows; the statistical
modules in all three tools can loop over rows. The current plan scopes to one row
(one pairwise comparison).

### 1.4 Reference Files

Only two reference files are required by this pipeline:

| Parameter | Description |
|-----------|-------------|
| `--gtf` * | Path to genome annotation in GTF format. **Must be the exact same GTF used in the nf-core/rnaseq run** that produced the input BAMs. Switching GTF sources (e.g., Ensembl for alignment, GENCODE for splicing) causes chromosome naming mismatches and silent failures. |
| `--transcript_fasta` | Transcript-level FASTA for IsoformSwitchAnalyzeR. **Optional when `--use_gffread true`** (strongly recommended). If supplied directly, it must originate from the same release and genome build as `--gtf`, because transcript ID version suffixes (e.g., `ENST00000001.1` in the FASTA vs `ENST00000001` in the GTF) cause silent import failures in ISAR. When `--use_gffread true`, this file is derived automatically from `--genome_fasta` + `--gtf`. |
| `--genome_fasta` | Genome-level FASTA. Required when `--use_gffread true`. Must match the genome sequence used in the upstream nf-core/rnaseq run. Used by gffread to extract transcript sequences that are guaranteed to be consistent with the provided GTF. |

### 1.5 Sequencing Data Quality Requirements

The following requirements apply to the input RNA-seq data regardless of input mode.
They should be verified before starting the pipeline and documented in the Methods section
of any publication.

#### Library preparation
- **Ribo-depletion is mandatory** for total RNA-Seq (RiboZero Gold, RNase H-based kits,
  or equivalent). Without ribo-depletion the vast majority of reads originate from rRNA
  and effective splicing coverage is insufficient. If poly-A selection was used instead,
  note that retained intron (RI) detection is inherently biased because pre-mRNA is
  largely excluded from poly-A-selected libraries.
- **Paired-end sequencing is strongly recommended** for splicing analysis. Single-end
  reads reduce junction-read sensitivity and are not supported by this pipeline.

#### Read length and sequencing depth
- Minimum recommended read length: **100 bp** (ideally 150 bp). Longer reads increase
  the number of junction-spanning reads and improve novel splice-site detection.
- Minimum recommended sequencing depth for alternative splicing analysis: **50 M
  paired-end reads per sample**. At depths below 30 M reads, many events will have
  insufficient junction coverage and will be filtered out before statistical testing,
  reducing sensitivity substantially.
- The `--read_length` parameter must reflect the **actual (post-trimming) read length**
  present in the BAM files, not the sequencing instrument nominal read length. The
  `validate_samplesheet.py` script uses **pysam** to extract the modal read length from
  the first aligned reads in each BAM and exits with a descriptive error if the
  declared value deviates by more than 10 bp. This is critical: both rMATS
  (`readLength` normalisation denominator) and MAJIQ (`readlen` config field) are
  sensitive to this value and will produce silently biased results if it is wrong.

#### RNA quality
- RNA Integrity Number (RIN) ≥ 7 is the minimum acceptable threshold. RIN ≥ 8 is
  preferred. Degraded RNA (RIN < 6) introduces systematic bias: degradation produces
  reads that map preferentially to the 3' end of transcripts and increases apparent
  intron retention because reads from incompletely processed pre-mRNA are more
  abundant relative to mature mRNA.
- Samples with RIN outliers (> 2 standard deviations below the group mean) should be
  excluded from the analysis or treated as a sensitivity check.

#### Chromosome naming consistency (GTF vs BAM)
This is the most common source of silent failure in splicing pipelines. The chromosome
names in the GTF annotation file **must exactly match** the chromosome names in the BAM
files.

- **Ensembl** style: `1, 2, 3, ..., X, Y, MT` (no `chr` prefix).
- **GENCODE / UCSC** style: `chr1, chr2, chr3, ..., chrX, chrY, chrM`.

Mixing sources (e.g., BAM aligned with an Ensembl genome but GTF from GENCODE) causes
rMATS to silently detect no events (no chromosome names match) and MAJIQ to build an
empty splice graph. The `validate_samplesheet.py` script should check that the `@SQ`
header lines of the first BAM in each group use the same naming convention as the
provided GTF. Specifically it should check whether the first sequence name in the
GTF starts with `chr` and whether the BAM `@SQ SN:` fields also start with `chr`,
and exit with an error if they differ.

#### Strandedness
Strandedness is a **global pipeline parameter** `--strandedness` (`unstranded`,
`forward`, or `reverse`). It applies to all samples — all samples from the same
experiment share the same library preparation protocol.

**Where to find the correct value**: check the nf-core/rnaseq MultiQC report →
`Strand specificity` section; or open any sample's
`star_salmon/<sample>/lib_format_counts.json` and read the `expected_format` field.

Do not guess or assume — an incorrect value causes rMATS to produce **silently wrong
PSI values for all genes on the negative strand**, which is roughly 50 % of genes.
MAJIQ will silently mis-assign reads as well.

Mapping from nf-core/rnaseq strandedness to splicing tool arguments:
| nf-core/rnaseq value | `--strandedness` | rMATS `--libType` | MAJIQ config `strandness` |
|----------------------|-----------------|-------------------|---------------------------|
| `unstranded` | `unstranded` | `fr-unstranded` | `none` |
| `forward` | `forward` | `fr-secondstrand` | `forward` |
| `reverse` | `reverse` | `fr-firststrand` | `reverse` |

---

## 2. Repository Structure (nf-core conventions)

```
alternative-splicing-nextflow/
├── assets/
│   ├── samplesheet.csv                # Example samplesheet
│   ├── comparisons.csv                # Example comparisons file
│   ├── multiqc_config.yml             # MultiQC configuration
│   └── schema_input.json              # JSON schema for samplesheet validation
├── bin/
│   ├── validate_samplesheet.py        # Entry-point samplesheet validator
│   ├── prepare_rmats_input.py         # Generates b1.txt/b2.txt or s1.txt/s2.txt
│   ├── prepare_majiq_config.py        # Generates majiq.conf ini file
│   └── render_report.R                # Master R Markdown rendering script
├── conf/
│   ├── base.config                    # Default resource profiles
│   ├── docker.config                  # Docker container definitions
│   ├── singularity.config             # Singularity container definitions
│   ├── test.config                    # Minimal test dataset parameters
│   ├── test_full.config               # Full-size test parameters
│   └── modules.config                 # Per-process resource overrides
├── docs/
│   ├── README.md                      # Main documentation
│   ├── usage.md                       # Usage instructions and all parameters
│   └── output.md                      # Description of all output files
├── modules/
│   ├── local/                         # Pipeline-specific processes (not in nf-core/modules)
│   │   ├── validate_input/
│   │   │   └── main.nf
│   │   ├── rmats_prep/
│   │   │   └── main.nf
│   │   ├── rmats_post/
│   │   │   └── main.nf
│   │   ├── majiq_build/
│   │   │   └── main.nf
│   │   ├── majiq_deltapsi/
│   │   │   └── main.nf
│   │   ├── majiq_voila_tsv/
│   │   │   └── main.nf
│   │   ├── gffread_transcriptome/
│   │   │   └── main.nf
│   │   ├── isar_import/
│   │   │   └── main.nf
│   │   ├── isar_switch_test/
│   │   │   └── main.nf
│   │   ├── isar_extract_orf/
│   │   │   └── main.nf
│   │   ├── isar_functional_annotation/
│   │   │   └── main.nf
│   │   ├── isar_switch_consequences/
│   │   │   └── main.nf
│   │   └── report/
│   │       └── main.nf
│   └── nf-core/                       # Installed from github.com/nf-core/modules.
│       │                              # Managed with `nf-core modules install/update`.
│       │                              # Versions tracked in modules.json at repo root.
│       └── multiqc/
│           └── main.nf
├── subworkflows/
│   └── local/
│       ├── input_check/
│       │   └── main.nf          # Validates samplesheet, emits per-sample channels
│       ├── rmats_analysis/
│       │   └── main.nf          # rMATS prep → post → filter
│       ├── majiq_analysis/
│       │   └── main.nf          # build → deltapsi → voila tsv
│       └── isoformswitchr_analysis/
│           └── main.nf          # All ISAR steps
├── workflows/
│   └── alternative_splicing.nf  # Top-level workflow, calls all subworkflows
├── main.nf                      # Pipeline entry point
├── nextflow.config              # Global config: profiles, params defaults
├── nextflow_schema.json         # JSON schema for all parameters
├── modules.json                 # nf-core tools registry: tracks installed module names,
│                                # their git SHA from nf-core/modules, and the branch used.
│                                # Updated automatically by `nf-core modules install/update`.
├── .nf-core.yml                 # nf-core lint configuration: declares pipeline type
│                                # ('pipeline'), org name, and which lint tests to skip.
├── CHANGELOG.md
├── CITATIONS.md                 # BibTeX / DOI citations for every tool used (required by
│                                # nf-core linting and expected by journal reviewers).
├── README.md
└── CHANGELOG.md
```

---

## 3. Pipeline Logic and Execution Flow

### 3.1 Top-Level DAG

The pipeline performs no pre-processing. Inputs arrive as pre-existing files from a
prior **nf-core/rnaseq** run.

```
     SAMPLESHEET VALIDATION
              │
    ┌─────────┴───────────────────────────────────────┐
    │                                           │
tuple(meta, bam, bai)               tuple(meta, salmon_dir)
    │                                           │
    │  BAM + GTF                  Salmon quant.sf + GTF
    │                                           │
┌───┴──────────┐                               │
│              │                               │
rMATS-turbo   MAJIQ V3           IsoformSwitchAnalyzeR
│              │                               │
└───┬──────────┘                               │
    └──────────────────────┬────────────────────┘
                           │
                     FINAL REPORT
                 (MultiQC + R Markdown)
```

### 3.2 Channel Strategy

All per-sample metadata (sample ID, condition, replicate) travel through the pipeline
attached as a `meta` map on every channel, following nf-core conventions:

```groovy
// meta map structure
meta = [id: 'ctrl_rep1', condition: 'control', replicate: 1]
```

Grouping into per-comparison channels (all samples of group1 vs. all samples of
group2) is performed using `groupTuple` keyed on a `comparison_id` that is constructed
in the `input_check` subworkflow by cross-referencing the samplesheet with the
comparisons file.

---

## 4. Module Specifications

> This pipeline contains **no pre-processing modules**. Trimming, STAR alignment, BAM
> indexing, and Salmon quantification are performed upstream by nf-core/rnaseq. The
> modules below cover only the three splicing analysis branches and the final report.

### 4.1 rMATS-turbo Branch

#### Overview

rMATS-turbo analyses five AS event types:
- SE: Skipped Exon
- A5SS: Alternative 5' Splice Site
- A3SS: Alternative 3' Splice Site
- MXE: Mutually Exclusive Exons
- RI: Retained Intron

The tool is run in **BAM mode** (rMATS `--b1`/`--b2` flags). The FASTQ mode of rMATS
(which runs STAR internally) is NOT used; alignment is handled explicitly by the
pipeline's own BAM input (not a separate alignment subworkflow).

**Novel splice-site detection** (`--novelSS`) is always enabled. This flag detects
AS events that involve at least one unannotated splice site and writes results to the
`fromGTF.novelSpliceSite.[AS_Event].txt` files in addition to the standard outputs.

The rMATS workflow is split into **prep** and **post** tasks to allow parallel
per-sample preprocessing when running on a cluster.

#### RMATS_PREP (process)

One process instance per BAM file.

- Container: `docker.io/xinglab/rmats-turbo:v4.3.0` (Bioconda package also available:
  `bioconda::rmats-turbo=4.3.0`)
- Command:

```bash
python rmats.py \
    --b1 ${bam_list_file} \
    --gtf ${gtf} \
    -t paired \
    --readLength ${params.read_length} \
    --variable-read-length \
    --allow-clipping \
    --libType ${rmats_lib_type} \
    --nthread ${task.cpus} \
    --od ${od} \
    --tmp ${tmp} \
    --task prep
```

Notes:
- `--b1` receives a single-line text file with one BAM path (one sample per prep run).
  The `bin/prepare_rmats_input.py` script generates these single-line files.
- `--variable-read-length` is always set to handle trimmed reads of variable length.
- `--allow-clipping` is always set to handle soft/hard-clipped reads from STAR.
- `--libType` maps from the `--strandedness` param: `unstranded → fr-unstranded`,
  `forward → fr-secondstrand`, `reverse → fr-firststrand`.
- Each prep run writes `*.rmats` files into its own `--tmp` directory.
- Output: `*.rmats` files collected per comparison.

#### RMATS_POST (process)

One process instance per pairwise comparison, after all per-sample prep jobs complete.

- Container: same as prep
- Command:

```bash
# Gather all .rmats files into a single tmp directory
python cp_with_prefix.py ${prefix} ${merged_tmp}/ ${prep_tmp_dir}/*.rmats

python rmats.py \
    --b1 ${b1_txt} \
    --b2 ${b2_txt} \
    --gtf ${gtf} \
    -t paired \
    --readLength ${params.read_length} \
    --variable-read-length \
    --allow-clipping \
    --libType ${rmats_lib_type} \
    --nthread ${task.cpus} \
    --tstat ${params.rmats_tstat_threads} \
    --cstat ${params.rmats_cstat} \
    --novelSS \
    --mil ${params.rmats_min_intron_length} \
    --mel ${params.rmats_max_exon_length} \
    --individual-counts \
    --od ${od} \
    --tmp ${merged_tmp} \
    --task post
```

Notes:
- `--b1` / `--b2` are text files listing all BAM paths for group1 and group2,
  respectively, separated by commas. Generated by `bin/prepare_rmats_input.py`.
- `--novelSS` is always enabled; the flag generates `fromGTF.novelSpliceSite.[AS].txt`
  output files.
- `--individual-counts` is always enabled to obtain per-replicate counts.
- `--cstat` default: `0.0001`; adjust via `--rmats_cstat` param.
- `--mil` default: `50` (min intron length for novel SS); adjust via
  `--rmats_min_intron_length`.
- `--mel` default: `500` (max exon length for novel SS); adjust via
  `--rmats_max_exon_length`.

#### rMATS Output Files Explained

All files are written to `results/rmats/<comparison_id>/`:

For each event type in `{SE, A5SS, A3SS, MXE, RI}`:

| File | Content |
|------|---------|
| `[AS].MATS.JC.txt` | Final results: junction-reads-only counting. Key columns: `GeneID`, `geneSymbol`, `chr`, `strand`, event coordinates, `IJC_SAMPLE_1`, `SJC_SAMPLE_1`, `IJC_SAMPLE_2`, `SJC_SAMPLE_2`, `IncLevel1`, `IncLevel2`, `IncLevelDifference`, `PValue`, `FDR`. |
| `[AS].MATS.JCEC.txt` | Same as JC but counting both junction and exon-body reads. |
| `fromGTF.[AS].txt` | All AS events detected from the GTF + RNA-seq data. |
| `fromGTF.novelJunction.[AS].txt` | Events from novel combinations of annotated splice sites (no unannotated edges). |
| `fromGTF.novelSpliceSite.[AS].txt` | Events involving at least one unannotated splice site (only when `--novelSS` is active). |
| `JC.raw.input.[AS].txt` | Raw junction read counts before statistical testing. |
| `JCEC.raw.input.[AS].txt` | Raw junction+exon counts before statistical testing. |
| `individualCounts.[AS].txt` | Per-replicate counts (enabled by `--individual-counts`). |
| `summary.txt` | Totals and significant-event counts at FDR ≤ 0.05. |

**Downstream filtering recommendation:**
- Primary filter: `FDR <= 0.05` AND `abs(IncLevelDifference) >= 0.1` (10% delta PSI).
- Use `MATS.JC.txt` as the primary result; `MATS.JCEC.txt` provides cross-validation.
- Novel events (from `fromGTF.novelSpliceSite`) should be interpreted with caution
  and reported separately.

---

### 4.2 MAJIQ V3 Branch

#### Overview

MAJIQ detects and quantifies Local Splicing Variations (LSVs). LSVs are general
representations of splicing that include binary junctions (analogous to SE, RI, etc.)
but also complex multiway splicing events that other tools may miss. MAJIQ V3
improves memory usage, runtime, and accuracy compared to V2.

MAJIQ does not assume that replicate counts are identical across the two groups
under comparison; the Bayesian framework estimates PSI distributions per replicate
and computes a deltaPSI posterior with an associated probability. This makes it
robust to outlier replicates and heterogeneous datasets.

The MAJIQ workflow has three sequential stages: build, deltapsi, and voila tsv.

#### MAJIQ_BUILD (process)

One process per comparison (all BAMs are processed together to build a shared splice
graph). Can be parallelised per sample if desired, but a single build over all samples
is preferred to define a unified junction set.

- Container: **custom Docker image** — MAJIQ is distributed under an academic licence
  and is not available on Biocontainers or Docker Hub. The repo must include a
  `containers/majiq/Dockerfile` that:
  1. Starts from `python:3.12-slim` (Python 3.12 is the supported version for V3)
  2. Installs `gcc`, `libhts-dev`, `zlib1g-dev` via apt
  3. Runs `pip install git+https://bitbucket.org/biociphers/majiq_academic.git`
  4. (Does **not** bake the licence file into the image — it is injected at runtime)

  The licence file is provided via the `MAJIQ_LICENSE_FILE` environment variable:
  ```groovy
  // In process MAJIQ_BUILD, MAJIQ_DELTAPSI, MAJIQ_VOILA_TSV:
  environment "MAJIQ_LICENSE_FILE=${params.majiq_license}"
  ```
  Users pass `--majiq_license /path/to/majiq_license*.key` at runtime.

- Configuration file generated by `bin/prepare_majiq_config.py`:

```ini
[info]
readlen = <read_length>
bamdirs = <comma-separated list of directories containing BAM files>
genome = <genome_build_string>  ; e.g., hg38 or mm10
strandness = <none|forward|reverse>

[experiments]
<sample_id_1> = <bam_filename_without_path>
<sample_id_2> = <bam_filename_without_path>
...
```

Notes on config:
- `readlen` must match `--read_length` param (same value used for rMATS).
- `strandness`: maps from `--strandedness` param: `unstranded → none`,
  `forward → forward`, `reverse → reverse`.
- `bamdirs` must be a comma-separated list of all directories that contain BAMs;
  since BAMs may live in different directories per sample (output from STAR), the
  prepare script will either symlink all BAMs into a shared directory or list each
  directory individually.

Command:

```bash
majiq build \
    ${gtf} \
    -c ${majiq_config} \
    -j ${task.cpus} \
    -o ${outdir} \
    --min-denovo-reads ${params.majiq_min_denovo_reads} \
    --min-intronic-cov ${params.majiq_min_intronic_cov}
```

Notes:
- `--min-denovo-reads` default: `2`; controls sensitivity for novel junction detection.
- `--min-intronic-cov` default: `3`; coverage threshold for retained intron detection.
- Outputs: one `.majiq` file per sample, `splicegraph.sql` (the unified splice graph),
  and `build.log`.

#### MAJIQ_DELTAPSI (process)

One process per pairwise comparison.

Command:

```bash
majiq deltapsi \
    -grp1 ${group1_majiq_files} \
    -grp2 ${group2_majiq_files} \
    --minreads ${params.majiq_min_reads} \
    --minnonzero ${params.majiq_min_nonzero} \
    -j ${task.cpus} \
    -n ${comparison_name_grp1} ${comparison_name_grp2} \
    -o ${outdir}
```

Notes:
- `-grp1` / `-grp2`: space-separated paths to `.majiq` files for each group (all
  replicates).
- `-n` provides human-readable labels for the output files.
- `--minreads` default: `10` (minimum total reads across all replicates to include an
  LSV).
- `--minnonzero` default: `3` (minimum number of replicates with at least one read).
- Outputs: `<comparison>.deltapsi.voila` file and `<comparison>.psi_grp1` /
  `<comparison>.psi_grp2` files (Bayesian PSI posteriors per group).

#### MAJIQ_VOILA_TSV (process)

Converts VOILA binary output to tab-separated text for downstream analysis.

Command:

```bash
voila tsv \
    ${splicegraph} \
    ${deltapsi_voila} \
    -f ${outdir} \
    --show-all \
    --threshold ${params.majiq_delta_psi_threshold} \
    --probability-threshold ${params.majiq_probability_threshold}
```

Notes:
- `--show-all`: outputs all LSVs, not just those passing the threshold; the threshold
  flags significant events without discarding data. Post-processing filters are applied
  in the report.
- `--threshold` default: `0.2` (20% minimum expected |deltaPSI|).
- `--probability-threshold` default: `0.95` (posterior probability that |deltaPSI| > threshold).
- Primary output: `<comparison_id>.tsv` with columns including `gene_name`, `lsv_id`,
  `junctions_coords`, `mean_dpsi_per_lsv_junction`, `probability_changing`,
  `probability_non_changing`, `group1_mean_psi`, `group2_mean_psi`.

#### MAJIQ Output Files Explained

All files written to `results/majiq/<comparison_id>/`:

| File | Content |
|------|---------|
| `splicegraph.sql` | Unified splice graph (SQLite) for all samples. Used by VOILA visualizer. |
| `*.majiq` | Per-sample Bayesian PSI posteriors (binary, MAJIQ internal format). |
| `*.deltapsi.voila` | Per-comparison deltaPSI posteriors (binary, VOILA input). |
| `<comparison>.tsv` | Tab-separated LSV results with PSI values and significance indicators. |

**Key interpretation notes:**
- MAJIQ reports LSVs, which may correspond to a single junction (simple binary events)
  or multiple junctions (complex events). The `lsv_id` encodes `gene:strand:position:type`.
- An LSV is considered "changing" when `probability_changing >= 0.95` AND
  `|mean_dpsi_per_lsv_junction| >= 0.2` (values from `--threshold` and
  `--probability-threshold` params).
- The absence of a classical p-value / FDR is by design; MAJIQ uses posterior
  probabilities from its Bayesian model instead.

---

### 4.3 IsoformSwitchAnalyzeR Branch

#### Overview

IsoformSwitchAnalyzeR (ISAR) identifies isoform switches — cases where two isoforms
of the same gene exchange ranking (i.e., one becomes more dominant while the other
becomes less dominant) between conditions. Starting from Salmon transcript-level
quantifications, ISAR performs a statistical test, extracts ORFs, and runs or imports
results from multiple external tools to annotate functional consequences.

The Salmon quantification used here must have been run with `--numBootstraps 100`
(covered in section 4.1) so that ISAR can use bootstrap uncertainty estimates rather
than point estimates, improving statistical accuracy.

ISAR is implemented as a series of R scripts orchestrated by Nextflow processes.
Each script is a self-contained R invocation that writes `.rds` objects as
intermediate results, which subsequent scripts read. This allows Nextflow to correctly
model the data flow.

#### GFFREAD_TRANSCRIPTOME (optional process)

When `params.use_gffread` is `true` (strongly recommended), this process runs before
`ISAR_IMPORT` and derives the transcript FASTA directly from the genome FASTA and the
GTF annotation. This guarantees 100% consistency between transcript IDs in the GTF and
in the FASTA, eliminating the silent import failures that arise when a separately
distributed Ensembl or GENCODE FASTA uses versioned transcript IDs
(e.g., `ENST00000001.1`) while the corresponding GTF uses unversioned IDs
(`ENST00000001`), or vice versa.

- Container: `biocontainers/gffread:0.12.7--hd03093a_1`
- Command:

```bash
gffread \
    ${gtf} \
    -g ${genome_fasta} \
    -w transcript_sequences.fa
```

- Inputs: `genome_fasta` (from `params.genome_fasta`), `gtf` (from `params.gtf`).
- Output: `transcript_sequences.fa` — passed directly to `ISAR_IMPORT` as the
  `isoformNtFasta` source, replacing the user-provided `--transcript_fasta` when
  `params.use_gffread` is `true`.

#### ISAR_IMPORT (process)

- Container: `bioconductor/bioconductor_docker:RELEASE_3_19` with ISAR and dependencies
  pre-installed. A custom Docker image is required:
  ```
  FROM bioconductor/bioconductor_docker:RELEASE_3_19
  RUN R -e "BiocManager::install('IsoformSwitchAnalyzeR')"
  RUN R -e "devtools::install_github('kvittingseerup/pfamAnalyzeR')"
  ```
  Image tag: `local/isar:2.10.0` (build and push to registry of choice).
  
- Inputs:
  - All `quant.sf` paths from Salmon (all samples, both conditions).
  - Samplesheet (used to construct the `colData` / design data frame).
  - GTF file.
  
- R script logic (`bin/isar_import.R`):

```r
library(IsoformSwitchAnalyzeR)
library(tximeta)

# Construct colData from samplesheet CSV
coldata <- read.csv(snakemake@input[["samplesheet"]])  # adapt for Nextflow param passing
coldata <- data.frame(
    sampleID  = coldata$sample,
    condition = coldata$condition
)

# Import Salmon quantifications using tximeta-aware import
switchAnalyzeRlist <- importRdata(
    isoformCountMatrix   = NULL,    # use importIsoformExpression below
    isoformRepExpression = NULL,
    designMatrix         = coldata,
    isoformExonAnnoation = snakemake@input[["gtf"]],  # GTF path
    isoformNtFasta       = snakemake@input[["transcript_fasta"]],
    showProgress         = FALSE
)
# Alternative: use importSalmonData() when Salmon output dirs are available
# which is the preferred path
switchAnalyzeRlist <- importSalmonData(
    sampleAnnotations    = coldata,  # data.frame with sampleID and condition
    salmonFileDir        = NULL,     # set via vector of quant.sf paths
    addIFmatrix          = TRUE,
    addBootstraps        = TRUE      # requires --numBootstraps in Salmon
)
saveRDS(switchAnalyzeRlist, file = "switchAnalyzeRlist_raw.rds")
```

- Output: `switchAnalyzeRlist_raw.rds`

#### ISAR_SWITCH_TEST (process)

Performs the statistical isoform switch test using DEXSeq or satuRn (default: satuRn
for improved accuracy with small sample sizes, as recommended in ISAR ≥ 2.x).

- R script logic:

```r
switchAnalyzeRlist <- readRDS("switchAnalyzeRlist_raw.rds")

# Filter isoforms before testing
switchAnalyzeRlist <- preFilter(
    switchAnalyzeRlist,
    geneExpressionCutoff = 1,       # min FPKM/TPM across samples
    isoformExpressionCutoff = 0,
    IFcutoff = 0.01,                # min isoform fraction
    removeSingleIsoformGenes = TRUE,
    reduceToSwitchingGenes = FALSE  # keep all for now; filter post-test
)

# Statistical test: satuRn (recommended for n < 10 replicates per group)
switchAnalyzeRlist <- isoformSwitchTestSatuRn(
    switchAnalyzeRlist,
    reduceToSwitchingGenes = TRUE,
    alpha                  = 0.05,
    dIFcutoff              = 0.1    # min 10% change in isoform fraction
)

saveRDS(switchAnalyzeRlist, "switchAnalyzeRlist_tested.rds")
```

- Output: `switchAnalyzeRlist_tested.rds`

#### ISAR_EXTRACT_ORF (process)

Extracts open reading frames from transcript sequences; required for functional
annotation steps.

```r
switchAnalyzeRlist <- readRDS("switchAnalyzeRlist_tested.rds")

switchAnalyzeRlist <- addORFfromGTF(
    switchAnalyzeRlist,
    pathToGTF = params$gtf,
    overwriteExistingORF = TRUE
)

# Extract ORFs from transcripts not annotated in GTF
switchAnalyzeRlist <- analyzeNovelIsoformORF(
    switchAnalyzeRlist,
    analysisAllIsoformsWithoutORF = TRUE,
    genomeObject = NULL   # genome BSgenome NOT required when transcript FASTA is provided
)

# Export FASTA sequences for external annotation tools
exportSequences(
    switchAnalyzeRlist,
    pathToOutput       = "./",
    writeTranscriptSequences = TRUE,
    writeORFSequences = TRUE,
    writePeptideSequences = TRUE
)

saveRDS(switchAnalyzeRlist, "switchAnalyzeRlist_orf.rds")
```

- Outputs: `switchAnalyzeRlist_orf.rds`, `isoformSwitchAnalyzeR_isoform.fasta`,
  `isoformSwitchAnalyzeR_ORF.fasta`, `isoformSwitchAnalyzeR_AA.fasta`
  (AA FASTA required for external tools).

#### ISAR_FUNCTIONAL_ANNOTATION (process)

This process runs external functional annotation tools on the exported amino acid FASTA
and imports the results back into the IsoformSwitchAnalyzeR object. The external tools
required are:

| Tool | Purpose | Container |
|------|---------|-----------|
| CPC2 | Coding potential prediction | `bioconda::cpc2` |
| PFAM (via Hmmer) | Protein domain annotation | `bioconda::hmmer` + Pfam-A.hmm |
| SignalP 6.0 | Signal peptide prediction | `brunopontecvo/signalp6` |
| DeepTMHMM | Transmembrane helix prediction | `biolib/deeptmhmm` |
| IUPred3 | Intrinsically disordered regions (IDR) | web API or local install |
| NetSurfP-3.0 | Secondary structure + solvent acc. | `biolib/netsurfp-3` |

Since several of these tools (SignalP, DeepTMHMM, NetSurfP) require BioLib or
institutional licences, they should be run as individual `exec` commands within a
single Nextflow process using a combined container image, OR as separate processes
with the results collected and passed together to `ISAR_SWITCH_CONSEQUENCES`.

The R import commands that follow each external run:

```r
# After CPC2
switchAnalyzeRlist <- analyzeCPC2(
    switchAnalyzeRlist,
    pathToCPC2resultFile = "cpc2_output.txt",
    removeNoncodinORFs   = TRUE
)

# After PFAM (Hmmer search against Pfam-A.hmm, parsed to PFAM format)
switchAnalyzeRlist <- analyzePFAM(
    switchAnalyzeRlist,
    pathToPFAMresultFile = "pfam_results.txt",
    showProgress = FALSE
)

# After SignalP 6.0
switchAnalyzeRlist <- analyzeSignalP(
    switchAnalyzeRlist,
    pathToSignalPresultFile = "signalp6_output.txt"
)

# After DeepTMHMM
switchAnalyzeRlist <- analyzeDeepTMHMM(
    switchAnalyzeRlist,
    pathToDeepTMHMMresultFile = "deeptmhmm_output.txt"
)

# After IUPred3
switchAnalyzeRlist <- analyzeIUPred2A(
    switchAnalyzeRlist,
    pathToIUPred2AresultFile = "iupred3_output.txt",
    extractorType            = "long",   # IDR predictions
    showProgress             = FALSE
)

# After NetSurfP-3.0
switchAnalyzeRlist <- analyzeNetSurfP2(
    switchAnalyzeRlist,
    pathToNetSurfP2resultFile = "netsurfp3_output.txt"
)

saveRDS(switchAnalyzeRlist, "switchAnalyzeRlist_annotated.rds")
```

- Output: `switchAnalyzeRlist_annotated.rds`

#### ISAR_SWITCH_CONSEQUENCES (process)

Analyses and summarises which functional features are gained or lost upon isoform
switching.

```r
switchAnalyzeRlist <- readRDS("switchAnalyzeRlist_annotated.rds")

# Compute functional consequences
switchAnalyzeRlist <- analyzeSwitchConsequences(
    switchAnalyzeRlist,
    consequencesToAnalyze = c(
        'intron_retention',
        'coding_potential',
        'ORF_seq_similarity',
        'NMD_status',
        'domains_identified',
        'domain_isotype',
        'IDR_identified',
        'IDR_type',
        'IDR_seq_similarity',
        'signal_peptide_identified',
        'signal_peptide_confidence',
        'topology_identified',
        'extracellular_domain'
    ),
    dIFcutoff  = 0.1,
    onlySigIsoforms = FALSE,
    showProgress    = FALSE
)

# Extract top switching genes
topSwitches <- extractTopSwitches(
    switchAnalyzeRlist,
    filterForConsequences = TRUE,
    n                     = Inf,     # return all, filter in report
    extractGenes          = FALSE,
    sortByQvals           = TRUE
)
write.csv(topSwitches, "top_isoform_switches.csv", row.names = FALSE)

# Extract consequence summary
consequenceSummary <- extractConsequenceSummary(
    switchAnalyzeRlist,
    includeCombined    = TRUE,
    consequencesToPlot = "all",
    asFractionTotal    = FALSE
)
write.csv(consequenceSummary, "consequence_summary.csv", row.names = FALSE)

# Generate switch plots for top N genes
switchPlotTopN(
    switchAnalyzeRlist,
    n = 25,
    filterForConsequences = TRUE,
    fileType = "pdf",
    pathToOutput = "switchplots/"
)

saveRDS(switchAnalyzeRlist, "switchAnalyzeRlist_final.rds")
```

- Outputs:
  - `switchAnalyzeRlist_final.rds` — final R object with all annotations and results
  - `top_isoform_switches.csv` — table of significant switches with consequences
  - `consequence_summary.csv` — genome-wide summary of consequence types
  - `switchplots/` — PDF plots per gene showing isoform structure, expression,
    and functional annotation tracks

#### ISAR Output Files Explained

All files written to `results/isoformswitchr/<comparison_id>/`:

| File | Content |
|------|---------|
| `switchAnalyzeRlist_final.rds` | Complete R object; load in R for further custom analysis. |
| `top_isoform_switches.csv` | Per-isoform listing: gene name, isoform IDs, dIF, q-value, consequence types. `dIF` = difference in isoform fraction (analogous to deltaPSI). |
| `consequence_summary.csv` | Counts of each consequence type (domain gain/loss, signal peptide gain/loss, IDR change, NMD sensitivity change, etc.) observed across all significant switches. |
| `switchplots/*.pdf` | Multi-panel plot per gene: transcript structure | protein domains | IDR | signal peptide | expression bar chart | isoform fraction bar chart. |
| `isoformSwitchAnalyzeR_AA.fasta` | Amino acid sequences for all tested isoforms (input to external annotation tools). |

---

## 5. Final Report

### 5.1 MultiQC

MultiQC is run at the end of the pipeline and collects logs from:
- FastQC and trimming QC (if available from nf-core/rnaseq MultiQC JSON)
- STAR alignment log (`Log.final.out` from BAM header, if available)
- rMATS `summary.txt`

Configuration file: `assets/multiqc_config.yml`
Output: `results/multiqc/multiqc_report.html`

### 5.2 R Markdown Report

A master R Markdown document (`bin/render_report.Rmd`) is rendered by a dedicated
Nextflow process and produces a single HTML report (with floating TOC).

The report has the following sections:

1. **Summary** — brief description of the experiment, number of samples per group,
   tools run, filter thresholds applied.
   
2. **Quality Control** — embedded MultiQC plots: mapping rates,
   duplication, insert size, gene body coverage.
   
3. **rMATS-turbo Results**
   - Total events detected per AS type (annotated, novel junction, novel splice site).
   - Volcano-style plot: IncLevelDifference (x-axis) vs. -log10(FDR) (y-axis), per AS type.
   - Bar chart: count of significant events (FDR ≤ 0.05 & |dPSI| ≥ 0.1) per event type.
   - Table: top 20 significant events per event type (sorted by FDR).
   - Venn diagram or UpSet plot: overlap of significant genes across event types.
   - Separate section for novel splice site events.

4. **MAJIQ V3 Results**
   - Total LSVs detected.
   - Distribution of |meandeltaPSI| for all detected LSVs.
   - Count of LSVs meeting the `P(|dPSI| > threshold`) ≥ 0.95 criterion.
   - Table: top 20 changing LSVs sorted by `probability_changing`.
   - Splice-type composition (inferred from LSV ID: source vs. target LSV).

5. **IsoformSwitchAnalyzeR Results**
   - Number of isoform switches identified (at dIF ≥ 0.1, q ≤ 0.05).
   - Consequence summary bar chart (domain, IDR, signal peptide, NMD, topology).
   - Table: top 25 switches with functional consequences.
   - Embed switch plots for the top 5 genes.

6. **Cross-Tool Overlap**
   - UpSet or Venn diagram: genes significant in rMATS AND MAJIQ AND ISAR.
   - Interpretation guidance: convergent hits from multiple tools are the highest
     confidence candidates.

7. **Methods** — auto-generated text listing tool versions, parameters, and reference.

8. **Session Info** — R session info for reproducibility.

Report rendering command (called from Nextflow `report` process):

```bash
Rscript -e "
  rmarkdown::render(
    input        = '${report_rmd}',
    output_format = rmarkdown::html_document(
      toc = TRUE, toc_float = TRUE, toc_depth = 3,
      code_folding = 'hide', theme = 'flatly'
    ),
    output_file  = '${output_html}',
    params       = list(
      rmats_dir   = '${rmats_dir}',
      majiq_dir   = '${majiq_dir}',
      isar_rds    = '${isar_rds}',
      comparison  = '${comparison_id}',
      fdr_cutoff  = ${params.report_fdr_cutoff},
      dpsi_cutoff = ${params.report_dpsi_cutoff}
    )
  )
"
```

---

## 6. Pipeline Parameters (nextflow_schema.json entries)

All parameters are declared in `nextflow_schema.json`. The following table covers every
parameter the pipeline exposes. Parameters with a `*` are required; others have defaults.

### Input / Output

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--input` * | string | — | Path to samplesheet CSV. |
| `--comparisons` * | string | — | Path to comparisons CSV. |
| `--outdir` | string | `./results` | Results directory. |

### Reference Files

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--gtf` * | string | — | GTF annotation file. Must be the same GTF used in the upstream nf-core/rnaseq run. |
| `--transcript_fasta` | string | — | Transcript FASTA for IsoformSwitchAnalyzeR. Required unless `--use_gffread true`. Must originate from the same release and build as `--gtf`. |
| `--genome_fasta` | string | — | Genome-level FASTA. Required when `--use_gffread true`. Must match the genome used in the upstream nf-core/rnaseq run. |
| `--use_gffread` | boolean | `false` | When `true`, run gffread to extract the transcript FASTA from `--genome_fasta` + `--gtf` before ISAR import. Strongly recommended to eliminate transcript ID version mismatches. Requires `--genome_fasta`. |

### Library Properties

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--strandedness` * | string | — | `unstranded`, `forward`, or `reverse`. Check nf-core/rnaseq MultiQC report (`Strand specificity` section) or `lib_format_counts.json` for the correct value. |
| `--read_length` * | integer | — | Read length (after trimming) used in the nf-core/rnaseq run. Required by rMATS and MAJIQ. |
| `--is_gencode` | boolean | `false` | Whether GTF is from GENCODE (affects ISAR transcript FASTA import). |

### rMATS parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--run_rmats` | boolean | `true` | Enable/disable rMATS branch. |
| `--rmats_cstat` | number | `0.0001` | Cutoff for null hypothesis test in rMATS. |
| `--rmats_min_intron_length` | integer | `50` | Min intron length for novel SS. |
| `--rmats_max_exon_length` | integer | `500` | Max exon length for novel SS. |
| `--rmats_tstat_threads` | integer | `6` | Threads for rMATS statistical model. |

### MAJIQ parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--run_majiq` | boolean | `true` | Enable/disable MAJIQ branch. |
| `--majiq_license` | path | — | **Required.** Path to the MAJIQ academic licence file (file whose name begins with `majiq_license`). Injected as `MAJIQ_LICENSE_FILE` inside the container. |
| `--majiq_min_denovo_reads` | integer | `2` | Min reads for de novo junction detection. |
| `--majiq_min_intronic_cov` | integer | `3` | Min intronic coverage for RI detection. |
| `--majiq_min_reads` | integer | `10` | Min reads per LSV for deltapsi. |
| `--majiq_min_nonzero` | integer | `3` | Min replicates with reads for deltapsi. |
| `--majiq_delta_psi_threshold` | number | `0.2` | deltaPSI threshold for significance. |
| `--majiq_probability_threshold` | number | `0.95` | Posterior probability threshold. |

### IsoformSwitchAnalyzeR parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--run_isar` | boolean | `true` | Enable/disable ISAR branch. |
| `--isar_alpha` | number | `0.05` | Adjusted p-value cutoff for switch test. |
| `--isar_dif_cutoff` | number | `0.1` | Min |dIF| for isoform switch. |
| `--isar_gene_expr_cutoff` | number | `1` | Min mean gene expression (TPM). |

### Report parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--report_fdr_cutoff` | number | `0.05` | FDR cutoff for rMATS event summary in report. |
| `--report_dpsi_cutoff` | number | `0.1` | |deltaPSI| cutoff for rMATS event summary in report. |

---

## 7. Configuration Files

### nextflow.config (structure)

```groovy
manifest {
    name            = 'alternative-splicing-nextflow'
    author          = '<author>'
    homePage        = 'https://github.com/<org>/alternative-splicing-nextflow'
    description     = 'Multi-tool alternative splicing Nextflow pipeline'
    mainScript      = 'main.nf'
    nextflowVersion = '!>=23.10.0'
    version         = '1.0.0'
}

params {
    // Input/output
    input        = null
    comparisons  = null
    outdir       = './results'
    // Reference files (must match what was used in the upstream nf-core/rnaseq run)
    // NOTE: pass via CLI or -params-file yaml, never via -c config files.
    // See: https://nf-co.re/docs/usage/configuration#custom-configuration-files
    gtf              = null
    genome_fasta     = null   // required when use_gffread = true
    transcript_fasta = null   // optional when use_gffread = true
    use_gffread      = false  // derive transcript FASTA from genome_fasta + gtf via gffread
    // Library properties (from the nf-core/rnaseq run)
    strandedness = null     // 'unstranded', 'forward', or 'reverse'
    read_length  = null     // read length after trimming
    is_gencode   = false
    // Tool switches
    run_rmats     = true
    run_majiq     = true
    run_isar      = true
    // rMATS defaults
    rmats_cstat               = 0.0001
    rmats_min_intron_length   = 50
    rmats_max_exon_length     = 500
    rmats_tstat_threads       = 6
    // MAJIQ (licence file required — obtain from majiq.biociphers.org)
    majiq_license             = null  // path to majiq_license*.key file
    // MAJIQ defaults
    majiq_min_denovo_reads    = 2
    majiq_min_intronic_cov    = 3
    majiq_min_reads           = 10
    majiq_min_nonzero         = 3
    majiq_delta_psi_threshold = 0.2
    majiq_probability_threshold = 0.95
    // ISAR defaults
    isar_alpha            = 0.05
    isar_dif_cutoff       = 0.1
    isar_gene_expr_cutoff = 1
    // Report defaults
    report_fdr_cutoff  = 0.05
    report_dpsi_cutoff = 0.1
}

profiles {
    docker {
        docker.enabled         = true
        docker.runOptions      = '-u $(id -u):$(id -g)'
        singularity.enabled    = false
    }
    singularity {
        singularity.enabled    = true
        singularity.autoMounts = true
        docker.enabled         = false
    }
    test {
        includeConfig 'conf/test.config'
    }
    test_full {
        includeConfig 'conf/test_full.config'
    }
}

includeConfig 'conf/base.config'
includeConfig 'conf/modules.config'
```

### conf/base.config

Defines default resource allocations and retry logic:

```groovy
process {
    cpus   = { check_max(2 * task.attempt, 'cpus') }
    memory = { check_max(8.GB * task.attempt, 'memory') }
    time   = { check_max(4.h * task.attempt, 'time') }

    errorStrategy = { task.exitStatus in [143, 137, 104, 134, 139] ? 'retry' : 'finish' }
    maxRetries    = 3
    maxErrors     = '-1'
}
```

### conf/modules.config

Per-process resource overrides. Suggested starting values:

```groovy
process {
    withName: 'STAR_ALIGN' {
        cpus   = 16
        memory = 64.GB
        time   = 8.h
    }
    withName: 'RMATS_PREP' {
        cpus   = 4
        memory = 16.GB
        time   = 4.h
    }
    withName: 'RMATS_POST' {
        cpus   = 8
        memory = 32.GB
        time   = 8.h
    }
    withName: 'MAJIQ_BUILD' {
        cpus   = 8
        memory = 32.GB
        time   = 6.h
    }
    withName: 'MAJIQ_DELTAPSI' {
        cpus   = 8
        memory = 16.GB
        time   = 4.h
    }
    withName: 'ISAR_.*' {
        cpus   = 4
        memory = 32.GB
        time   = 4.h
    }
    withName: 'SALMON_QUANT' {
        cpus   = 8
        memory = 16.GB
        time   = 4.h
    }
}
```

---

## 8. Container Strategy

Each tool should have its own pinned container image. The pipeline uses container
scoping per process (`process.container` in modules) to allow mixing Docker and
Singularity across processes.

| Process(es) | Image |
|-------------|-------|
| FastQC | `biocontainers/fastqc:v0.12.1_cv4` | upstream (nf-core/rnaseq) |
| Trim Galore | `biocontainers/trim-galore:0.6.10--hdfd78af_0` | upstream (nf-core/rnaseq) |
| STAR | `biocontainers/star:2.7.10a--h9ee0642_0` |
| Samtools | `biocontainers/samtools:1.19.2--h50ea8bc_0` |
| Salmon | `biocontainers/salmon:1.10.1--h0a6ea8b_1` |
| rMATS-turbo | `xinglab/rmats-turbo:v4.3.0` |
| MAJIQ V3 | **custom image** — `containers/majiq/Dockerfile` in this repo. Install via `pip install git+https://bitbucket.org/biociphers/majiq_academic.git` (Python 3.12 + GCC + HTSlib + zlib). Licence file injected at runtime via `MAJIQ_LICENSE_FILE` env var. |
| gffread | `biocontainers/gffread:0.12.7--hd03093a_1` |
| IsoformSwitchAnalyzeR | custom image based on `bioconductor/bioconductor_docker:RELEASE_3_19` |
| PFAM + Hmmer | `biocontainers/hmmer:3.3.2--hdbdd923_4` |
| SignalP 6 | institutional licence required; pull from BioLib if available |
| DeepTMHMM | `biolib/deeptmhmm:latest` or BioLib API |
| IUPred3 | `bioconda::iupred` or local install |
| NetSurfP-3.0 | `biolib/netsurfp-3:latest` or BioLib API |
| R Markdown report | `rocker/verse:4.3.2` (includes rmarkdown, ggplot2, tidyverse) |
| MultiQC | `biocontainers/multiqc:1.21--pyhdfd78af_0` |

**Note on SignalP / DeepTMHMM / NetSurfP licencing:** These three tools are academic
licences. They can be used offline by downloading the model weights and running
locally, or via BioLib cloud calls. The pipeline should accept pre-computed result
files via parameters (`--signalp_results`, `--deeptmhmm_results`, `--netsurfp_results`)
as an alternative to running these tools internally, so that the pipeline does not
block on licencing availability.

---

## 9. Output Directory Structure

```
results/
├── qc/
│   └── multiqc/
│       └── multiqc_report.html
├── rmats/
│   └── <comparison_id>/                 # e.g., control_vs_treatment
│       ├── SE.MATS.JC.txt
│       ├── SE.MATS.JCEC.txt
│       ├── fromGTF.SE.txt
│       ├── fromGTF.novelJunction.SE.txt
│       ├── fromGTF.novelSpliceSite.SE.txt
│       ├── JC.raw.input.SE.txt
│       ├── JCEC.raw.input.SE.txt
│       ├── individualCounts.SE.txt
│       ├── [same set for A5SS, A3SS, MXE, RI]
│       └── summary.txt
├── majiq/
│   └── <comparison_id>/
│       ├── splicegraph.sql
│       ├── build.log
│       ├── <sample_id>.majiq          # one per sample
│       ├── <comparison_id>.deltapsi.voila
│       └── <comparison_id>.tsv
├── isoformswitchr/
│   └── <comparison_id>/
│       ├── switchAnalyzeRlist_final.rds
│       ├── top_isoform_switches.csv
│       ├── consequence_summary.csv
│       ├── isoformSwitchAnalyzeR_AA.fasta
│       └── switchplots/
│           └── <gene_id>_switchplot.pdf
└── report/
    └── <comparison_id>_splicing_report.html
```

---

## 10. Implementation Roadmap

The following steps should be followed in order when implementing the pipeline.

### Step 1. Repository Initialisation
1. Create the directory tree exactly as described in Section 2.
2. Initialise a `CHANGELOG.md` and `README.md` with placeholders.
3. Create `nextflow.config` following the skeleton in Section 7.
4. Create `nextflow_schema.json` with all parameters from Section 6, following the
   nf-core JSON Schema specification (see https://nf-co.re/tools#pipeline-schema).
5. Create `.nf-core.yml` with minimal content:
   ```yaml
   nf_core_version: 2.13.0
   org: <your-org>
   template:
     name: alternative-splicing-nextflow
     version: 1.0.0
   lint:
     nextflow_config:
       - manifest.name
       - manifest.homePage
   ```
6. Create an empty `modules.json` with the correct skeleton so that `nf-core modules`
   can populate it:
   ```json
   {
     "name": "alternative-splicing-nextflow",
     "homePage": "https://github.com/<org>/alternative-splicing-nextflow",
     "repos": {
       "https://github.com/nf-core/modules.git": {
         "modules": {},
         "subworkflows": {}
       }
     }
   }
   ```
7. Add `.gitignore` (standard nf-core template: exclude `.nextflow/`, `work/`,
   `results/`, `*.log`).
8. Create `CITATIONS.md` with citation stubs for every tool listed in Section 12.

### Step 2. Input Validation Subworkflow
1. Implement `bin/validate_samplesheet.py` — reads the CSV, validates column names,
   checks that `bam`, `bai`, and `salmon_dir` paths exist for every row, checks that
   exactly two conditions are present, and exits with a human-readable error if
   validation fails. Also validates that `bam` chromosome names match the `--gtf`
   chromosome convention (`chr`-prefixed vs. bare) to catch silent mismatches early.
   Additionally, uses **pysam** to sample the first aligned reads from each BAM and
   extract the modal read length. This value is compared against the user-supplied
   `--read_length` parameter; if the difference exceeds 10 bp the script exits with a
   descriptive error identifying the discrepant sample(s). This prevents rMATS from
   silently producing biased PSI estimates and prevents MAJIQ from using an incorrect
   `readlen` value in its configuration file.
2. Implement the `input_check` subworkflow (`subworkflows/local/input_check/main.nf`)
   that calls the validation script and emits two channels:
   - `ch_samples_bam` — tuple(meta, bam, bai) for all samples
   - `ch_samples_salmon` — tuple(meta, salmon_dir) for all samples
3. Cross-reference with comparisons CSV to attach `comparison_id` and `group` to each
   sample in `meta`.

### Step 3. rMATS Subworkflow
1. Implement `bin/prepare_rmats_input.py`:
   - Accepts: list of BAM paths and their group membership.
   - Outputs: `b1.txt` and `b2.txt` files (comma-separated BAM lists).
2. Implement `modules/local/rmats_prep/main.nf` — runs one prep per sample BAM.
3. Implement `modules/local/rmats_post/main.nf` — runs one post per comparison,
   collecting all prep `.rmats` files via `cp_with_prefix.py`.
4. Implement `subworkflows/local/rmats_analysis/main.nf`:
   - Use `groupTuple` on `comparison_id` to collect all BAMs per comparison.
   - Fan out RMATS_PREP across all BAMs.
   - After all prep jobs complete, run RMATS_POST.
   - Emit: `ch_rmats_results` — path to output directory per comparison.

### Step 4. MAJIQ Subworkflow
1. Implement `bin/prepare_majiq_config.py`:
   - Accepts: list of (sample_id, bam_path) tuples, read_length, strandedness,
     genome build string.
   - Outputs: `majiq.conf` ini file.
2. Implement `modules/local/majiq_build/main.nf`.
3. Implement `modules/local/majiq_deltapsi/main.nf`.
4. Implement `modules/local/majiq_voila_tsv/main.nf`.
5. Implement `subworkflows/local/majiq_analysis/main.nf` wiring the three steps.

### Step 5. IsoformSwitchAnalyzeR Subworkflow
1. Implement `modules/local/gffread_transcriptome/main.nf` (container:
   `biocontainers/gffread:0.12.7--hd03093a_1`). When `params.use_gffread` is `true`,
   this process runs before `ISAR_IMPORT` and produces the transcript FASTA from the
   genome FASTA and GTF. Wire the conditional in
   `subworkflows/local/isoformswitchr_analysis/main.nf`: if `params.use_gffread` run
   `GFFREAD_TRANSCRIPTOME` and pass its output FASTA to `ISAR_IMPORT`; otherwise pass
   `file(params.transcript_fasta)` directly.
2. Create the custom Docker image for ISAR (see Section 8). Write `Dockerfile` under
   `containers/isar/Dockerfile` and build it.
3. Implement each R script under `bin/`:
   - `isar_import.R`
   - `isar_switch_test.R`
   - `isar_extract_orf.R`
   - `isar_functional_annotation.R`
   - `isar_switch_consequences.R`
4. Implement corresponding modules under `modules/local/`.
5. Implement `subworkflows/local/isoformswitchr_analysis/main.nf` chaining all steps
   (including the conditional `GFFREAD_TRANSCRIPTOME` pre-step from item 1).
6. For external tools (SignalP, DeepTMHMM, etc.) that require licences, implement a
   conditional: if `--signalp_results` is provided, skip the SignalP run and import
   the pre-computed file; otherwise attempt to run the tool via BioLib or local binary.

### Step 6. Report
1. Write `bin/render_report.Rmd` — full R Markdown document with all sections from
   Section 5.2. Each section reads from the results directories passed as `params`.
2. Implement `modules/local/report/main.nf` — calls `Rscript -e rmarkdown::render(...)`.
3. Implement `modules/nf-core/multiqc/main.nf` — collects all QC channel outputs.

### Step 7. Top-Level Workflow
1. Implement `workflows/alternative_splicing.nf`:
   - Call `input_check` → emits `ch_samples_bam` and `ch_samples_salmon`.
   - Conditional on `params.run_rmats`: call `rmats_analysis` with `ch_samples_bam`.
   - Conditional on `params.run_majiq`: call `majiq_analysis` with `ch_samples_bam`.
   - Conditional on `params.run_isar`: call `isoformswitchr_analysis` with `ch_samples_salmon`.
   - Collect all MultiQC inputs.
   - Call `multiqc`.
   - Call `report`.
2. Implement top-level `main.nf` that includes `workflows/alternative_splicing.nf`.

### Step 8. Test Profile
1. Obtain or generate minimal test data (small BAM files from a public dataset, e.g.,
   a subset of a Snaptron or Recount3 study with known AS events, covering ≥ 2 replicates
   per group). These are the same BAMs that nf-core/rnaseq would produce; download them
   directly rather than re-running alignment.
2. Write `conf/test.config` pointing to these files, with all parameters set.
3. Run with `nextflow run main.nf -r 1.0.0 -profile test,docker` and verify the pipeline
   completes without errors.

   > **Reproducibility tip**: always pin the pipeline version with `-r <tag>` when
   > running on real data (`-r 1.0.0`, `-r 1.2.0`, etc.).

### Step 9. Documentation and Finalisation
1. Fill in `docs/usage.md`: full parameter descriptions, samplesheet format, and
   example run command. Include a note on `NXF_OPTS`:
   ```bash
   # Prevent out-of-memory errors in the Nextflow JVM (add to ~/.bashrc)
   export NXF_OPTS='-Xms1g -Xmx4g'
   ```
   Recommended way to pass parameters — via `-params-file yaml`:
   ```bash
   cat > params.yaml <<'EOF'
   input: /path/to/samplesheet.csv
   comparisons: /path/to/comparisons.csv
   outdir: ./results
   gtf: /ref/Homo_sapiens.GRCh38.110.gtf
   transcript_fasta: /ref/Homo_sapiens.GRCh38.110.transcripts.fa
   strandedness: reverse
   read_length: 150
   EOF

   nextflow run main.nf -r 1.0.0 -profile singularity -params-file params.yaml
   ```
2. Fill in `docs/output.md`: describe every output file (this document can serve as
   the source).
3. Write `CHANGELOG.md` entry for v1.0.0.
4. Complete `CITATIONS.md` with final DOIs and BibTeX entries for all tools.
5. Run `nf-core lint` and resolve all warnings before tagging v1.0.0:
   ```bash
   nf-core lint .
   ```
   Common lint checks: `nextflow_schema.json` completeness, `modules.json` consistency
   with installed modules, presence of `CITATIONS.md`, `test` profile availability,
   container pinning in all modules.
6. Tag the release and push: `git tag v1.0.0 && git push --tags`.

---

## 11. Key Implementation Notes and Pitfalls

### Data Quality and Experimental Design
- **Minimum 3 replicates per group is a hard biological requirement**, not a software
  constraint. With n = 2, the Bayesian model in rMATS and the posterior estimation in
  MAJIQ do not have sufficient information to distinguish genuine signal from sampling
  noise; false discovery rates will be unreliable.
- **Inter-sample RNA quality variability** (RIN differences within a group) can
  introduce spurious splicing signals. Before running the pipeline, inspect the gene body
  coverage plot produced by RSeQC or the nf-core/rnaseq MultiQC report: a strong 3' bias
  in a subset of samples indicates degradation and that sample may need to be excluded.
- **Samples with divergent coverage profiles** (e.g., a clear outlier in PCA of PSI
  values) should be investigated before accepting results. rMATS per-replicate PSI values
  (`IncLevel1` / `IncLevel2` columns, comma-separated) make it easy to spot outlier
  replicates for significant events.
- **Sequencing depth below 30 M PE reads** markedly reduces sensitivity. Events in
  lowly expressed genes will lack sufficient junction reads to pass internal coverage
  thresholds and will not be quantified, leading to an underestimate of the total number
  of differential AS events.

### Reference Genome and Annotation
- **GTF and genome FASTA must originate from the same source and release.** Mixing
  Ensembl and GENCODE produces chromosome naming mismatches (`1` vs `chr1`) that cause
  rMATS to find zero events silently and MAJIQ to build an empty splice graph.
- **The STAR index must be built with the same STAR version used for alignment.** STAR
  indices are not backward-compatible across major or minor versions. If the index was
  built with STAR 2.7.10a and alignment is attempted with STAR 2.7.11, the run will
  abort with a version mismatch error. Pin the STAR version in both `star_genomegenerate`
  and `star_align` processes to the same container image.
- **Chromosome naming must be verified programmatically** in `validate_samplesheet.py`:
  read the `@SQ` header lines from the provided BAM and compare `SN:` field prefixes
  against the first non-comment line of the GTF. Exit with a clear error message if they
  differ.

### rMATS
- The `--readLength` parameter must reflect the **actual (post-trimming) modal read
  length** present in the BAM files, not the sequencing instrument nominal read length.
  The pipeline validates this automatically: `validate_samplesheet.py` uses pysam to
  extract the observed read length from each BAM at startup and exits with an error if
  the declared `--read_length` deviates by more than 10 bp. If reads were trimmed to
  variable lengths in the nf-core/rnaseq run (the default with quality-based trimming),
  always pass `--variable-read-length` alongside `--readLength`. The `--readLength`
  value is used to compute `IncFormLen` and `SkipFormLen` normalisation denominators;
  an incorrect value shifts all PSI estimates systematically without any warning from
  rMATS itself.
- When running prep in parallel across many samples, `*.rmats` filenames may collide.
  Always use `cp_with_prefix.py` to merge temporary directories before running post.
- Always verify that the BAM paths given in `--b1`/`--b2` for the post step exactly
  match the full paths used in prep; relative paths will cause silent lookup failures.
- `fromGTF.novelSpliceSite.[AS].txt` files are only populated when `--novelSS` is
  set; if the file is empty, there are no novel splice sites in that event type for
  that comparison — this is an expected outcome for well-annotated genomes.
- **The `--tmp` directory can grow to 50–100 GB** for a mammalian genome analysis with
  8–10 samples at 60–80 M reads per sample. Provision temporary storage accordingly.
  On SLURM clusters, use local scratch (`/lscratch/$SLURM_JOB_ID`) to avoid filling
  shared filesystems.
- **Filter rMATS results with both FDR and |dPSI|.** The FDR alone is insufficient:
  with 50 000 detected SE events at FDR < 0.05, roughly 2 500 are expected false
  positives. Applying `FDR ≤ 0.05 AND |IncLevelDifference| ≥ 0.1` reduces this
  substantially. Additionally, requiring a minimum of 10 junction reads per event per
  replicate (using the `IJC` + `SJC` columns) removes low-coverage noise that the
  statistical model cannot reliably resolve.
- **Events in lowly expressed genes are less reliable.** Cross-reference significant
  rMATS events against DESeq2 normalized counts (if a parallel differential expression
  analysis was run) to confirm that the host gene has adequate expression in both groups.
- **`--cstat` is not a p-value cutoff.** The default value of `0.0001` means the null
  hypothesis is tested for any ΔΨ > 0.01 %. Run with the default and apply your own
  |dPSI| threshold at the post-processing stage in R for full flexibility.
- **Visual validation in IGV or UCSC Genome Browser** is strongly recommended for the
  top 10–20 events before reporting them. Load the BAMs directly and inspect the
  junction-spanning reads manually. This catches alignment artefacts (multimapping,
  repeat regions), confirms the correct interpretation of the event (inclusion vs.
  skipping direction), and provides publication-quality sashimi plot material.

### MAJIQ
- MAJIQ requires BAM files to be coordinate-sorted and indexed (`.bai`).
- The `bamdirs` entry in the config must be the **directory** path, not the BAM file
  path. All BAM files listed under `[experiments]` must be accessible in at least one
  of the listed `bamdirs`.
- For stranded libraries, `strandness` in the MAJIQ config must match the library
  type. Incorrect strandedness is a common source of reduced LSV detection.
- MAJIQ V3 uses GFF3 or GTF annotation. Confirm that the GTF parser is invoked
  correctly (MAJIQ accepts GTF directly but the flag may differ from V2; verify with
  the V3 documentation).

### IsoformSwitchAnalyzeR
- Salmon must be run with `--numBootstraps 100` to enable the bootstrap-aware
  uncertainty model in ISAR. Without bootstraps, only point estimates are used.
- `importSalmonData()` requires the sampleAnnotation `condition` column to exactly
  match between the two groups specified in the comparisons CSV.
- The `preFilter()` step significantly affects the number of switches detected.
  Start with default thresholds and report the number of isoforms and genes retained
  after filtering in the Methods section.
- External annotation tools (SignalP, DeepTMHMM) output formats must match the
  expected ISAR import format. ISAR's import functions parse specific column layouts;
  do not reformat these files.
- `analyzeIUPred2A()` in ISAR accepts IUPred3 output if the output format is
  compatible; verify the column names from the installed IUPred3 version match what
  ISAR expects.

### General Nextflow / nf-core
- All processes must declare `tag "${meta.id}"` for readable log output.
- All processes must declare `label` matching entries in `conf/modules.config`.
- Publish directories must be declared with `mode: 'copy'` for all final output
  files and `mode: 'symlink'` for intermediate files used only within the pipeline.
- Use `workflow.onComplete` in `main.nf` to log a summary with runtime, exit status,
  and `--outdir` path.
- Parameter validation must occur early (before any process runs) using the
  `validateParameters()` call from the nf-validation plugin.

---

## 12. External Tool Version Pinning

> Pre-processing tools (FastQC, Trim Galore, STAR, Samtools, Salmon) are not part of
> this pipeline. They are versioned and pinned within the **nf-core/rnaseq** upstream
> run. The versions below cover only tools installed in this repository.

| Tool | Version | Source |
|------|---------|--------|
| Nextflow | ≥ 23.10.0 | nextflow.io |
| nf-validation plugin | 1.1.3 | nf-core |
| FastQC | 0.12.1 | upstream: nf-core/rnaseq |
| Trim Galore | 0.6.10 | upstream: nf-core/rnaseq |
| STAR | 2.7.10a | Bioconda |
| Samtools | 1.19.2 | Bioconda |
| Salmon | 1.10.1 | Bioconda |
| rMATS-turbo | 4.3.0 | Docker Hub (xinglab) / Bioconda |
| MAJIQ | 3.0 | biociphers.org / Bioconda |
| IsoformSwitchAnalyzeR | 2.10.0 | Bioconductor 3.22 |
| pfamAnalyzeR | latest | GitHub (kvittingseerup) |
| Hmmer | 3.3.2 | Bioconda |
| Pfam-A.hmm | 37.0 | EMBL-EBI |
| MultiQC | 1.21 | Bioconda |
| R | 4.3.x | Rocker / Bioconductor |
| satuRn | ≥ 1.7.0 | Bioconductor |
| DEXSeq | current | Bioconductor |

---

## 13. Update Resilience — Why This Pipeline Stays Current

This section explains precisely why the architecture is stable and what the actual
maintenance surface looks like over time.

### The key distinction: nf-core/modules vs nf-core/rnaseq

This pipeline installs only **one** module from `github.com/nf-core/modules`: MultiQC.
All other nf-core/modules wrappers (FastQC, STAR, Salmon, Samtools, etc.) are used
exclusively inside nf-core/rnaseq, which is a completely separate monolith.

### What does actually require maintenance

| Event | Action required | Effort |
|-------|----------------|--------|
| nf-core/modules updates a module (e.g., MultiQC) | Optionally run `nf-core modules update multiqc`. Review changelog. | < 1 hour |
| A tool releases a security fix or major bug fix | Update the container tag in the relevant `modules/nf-core/<tool>/main.nf` and in the version table in Section 12. | < 1 hour |
| rMATS-turbo releases a new version | Update container tag in `modules/local/rmats_prep/main.nf` and `rmats_post/main.nf`. Test with `--profile test`. | 1–4 hours |
| MAJIQ releases V3.x → V4 with breaking CLI changes | Update `modules/local/majiq_build/main.nf`, `majiq_deltapsi/main.nf`, `majiq_voila_tsv/main.nf`, and `bin/prepare_majiq_config.py`. | 4–8 hours |
| IsoformSwitchAnalyzeR releases a breaking API change | Update R scripts in `bin/isar_*.R`. Run `--profile test` to catch failures. | 4–8 hours |
| Nextflow itself releases a new DSL change | Update `manifest.nextflowVersion` in `nextflow.config`. Run `nf-core lint`. | 1–2 hours |

### Good practices already followed

- All containers are pinned to exact version tags (no `latest`).
- `modules.json` tracks the exact git SHA of the `multiqc` nf-core module installed.
- `nf-core lint` is run before every release to catch API drift.
- The `test` profile provides a fast regression check after any update.
- Input/output contracts (samplesheet columns, output directory layout) are versioned
  and described in `nextflow_schema.json` and `docs/output.md`.
- The three analytical branches (rMATS, MAJIQ, ISAR) are fully independent subworkflows.
  Each can be toggled with `--run_rmats false` etc., so a breaking update in one tool
  does not block the use of the other two.

---

*End of pipeline planning document.*
