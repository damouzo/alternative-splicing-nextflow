# Output Files Documentation

This document describes all output files produced by the alternative-splicing-nextflow pipeline.

Outputs are organised so that the folder that gets shipped is physically separate
from the audit evidence (see `assets/results_schema.yaml`):

- **Shippable** — `results/deliverables/` (`README.md`, `run_info/`, `contrasts/` and
  `cross_contrast/`). This is what the recipient opens; it follows the results contract
  (`results_contract_version`) and is QA-checked by `validate_results_contract.py` at the
  end of every run. All paths *inside* it are relative to `deliverables/`, so the folder
  can be renamed, moved and zipped without breaking anything.
- **Audit-only** — `results/raw/` and `deliverables/run_info/pipeline_info/`. Native
  per-tool output and the Nextflow execution artefacts. Not needed to interpret the
  results; not part of the shipment.

Shipping is therefore: `mv deliverables <name> && zip -r <name>.zip <name>`.

The pipeline does not migrate or delete outputs from a previous run/layout, and
`publish_dir_mode=copy` never removes old files. Start from a clean `--outdir`
(or keep previous runs in a separate directory) so stale tables cannot end up in
`deliverables/`.

`--publish_level` controls how much of the audit-only `raw/` layer is written:

| Level      | `raw/` native tool dirs                | rMATS JCEC / individual counts |
|------------|----------------------------------------|--------------------------------|
| `core`     | not published                          | not published |
| `core_raw` | published (default)                    | only with `--publish_rmats_jcec` / `--publish_rmats_individual_counts` |
| `full`     | published                             | always published              |

Sashimi final PDFs are always a shippable deliverable (they live under
`deliverables/contrasts/<id>/plots/sashimi/`), independent of `publish_level`.

`--publish_raw=false` behaves like `core` for the `raw/` layer.

`publish_dir_mode` defaults to `copy`, so the shipped files are real files and not
symlinks into `work/`.

## Directory Structure

```
results/
├── deliverables/                    # SHIPPABLE: rename + zip to send
│   ├── README.md                    # what each folder is, how to read the TSVs
│   ├── run_info/                    # what was run and how to read the tables
│   │   ├── sample_index.tsv         # sample x comparison x group (no cluster paths)
│   │   ├── software_versions.yml    # single consolidated tool-version record
│   │   ├── qa_report.txt            # results-contract validation output
│   │   └── pipeline_info/           # Nextflow report, trace, timeline, DAG (audit-only)
│   ├── contrasts/
│   │   └── <comparison_id>/
│   │       ├── <comparison_id>_splicing_report.html
│   │       ├── tables/
│   │       │   ├── <c>.rmats.{master,significant,summary}.tsv
│   │       │   ├── <c>.majiq.{master,significant,summary}.tsv      (when run_majiq)
│   │       │   ├── <c>.isar.{master,significant,summary}.tsv       (when run_isar)
│   │       │   ├── <c>.leafcutter.{master,significant,summary}.tsv (when run_leafcutter)
│   │       │   ├── <c>.pegasas.{master,significant,summary}.tsv    (when run_pegasas)
│   │       │   ├── <c>.cross_tool.master.tsv                       (when build_cross_tool_master)
│   │       │   └── <c>.cross_tool.gene_summary.tsv
│   │       └── plots/
│   │           └── sashimi/
│   │               ├── sashimi_index.tsv                           (when run_sashimi)
│   │               └── <EVENT_TYPE>/*.pdf                          final sashimi PDFs
│   └── cross_contrast/
│       └── pegasas/                 # cross_contrast_summary.tsv, heatmap, UpSet
├── run_info/                        # run metadata, NOT shipped
│   ├── run_manifest.yaml            # pipeline, schema_version, params, tools, known issues
│   └── contrast_manifests/
│       └── <comparison_id>_contrast_manifest.yaml
└── raw/                             # audit-only, NOT shipped
    ├── rmats/<comparison_id>/…
    ├── majiq/<comparison_id>/…
    ├── isar/<comparison_id>/…
    ├── leafcutter/<comparison_id>/…
    ├── sashimi/<comparison_id>/sashimi_out/…   # intermediates (final PDFs are a deliverable)
    ├── pegasas/<comparison_id>/…    # + _shared/ pathway scores (cohort-level)
    └── _internal/sample_paths.tsv   # absolute input paths for audit
```

Manifests (`run_manifest.yaml`, per-contrast manifests) are run metadata and live in
`results/run_info/`, outside the shipped `deliverables/` tree.

## Shippable layer (`deliverables/`)

### run_info/

| File | Description |
|------|-------------|
| `sample_index.tsv` | `sample_id, condition, replicate, comparison_id, group` — deliberately free of absolute paths |
| `software_versions.yml` | every per-process `versions.yml` merged into one file |
| `qa_report.txt` | output of `validate_results_contract.py` (errors + warnings) |
| `pipeline_info/` | Nextflow `execution_report.html`, `execution_trace.txt`, `execution_timeline.html`, `pipeline_dag.svg` (audit-only) |

`run_manifest.yaml` and the per-contrast manifests are in `results/run_info/` (outside
`deliverables/`), together with the embedded `params_criticos`, `tools_enabled`/
`tools_disabled`, `known_issues` and per-contrast manifests.

The absolute BAM/BAI/Salmon paths live only in `raw/_internal/sample_paths.tsv`, so they
do not leak into the shippable layer.

### Master table standard columns

Every `<tool>.master.tsv` starts with:

`comparison_id, tool, feature_type, feature_id, gene_id, gene_symbol, group1_name,
group2_name, effect_size_direction, effect_size, effect_size_type, pvalue, padj,
padj_method, is_significant, significance_rule, source_file`

followed by the tool-native columns (never dropped). `effect_size` is always
`group2 - group1` and `effect_size_direction` is `group2_minus_group1`. `padj` is always
accompanied by `padj_method` (`BH`, `rmats_cstat`, `empirical_FDR`, `BH_event_wise` or
`none`) so nobody filters on an ambiguous empty `fdr`. rMATS tables add
`event_type/event_id`, `fdr_floor_flag` (FDR==0 rows are the numeric floor of
`--cstat`), `event_class`, `reads_ok` and `is_igv_supported`.

rMATS location columns (contract 2.1.0): `event_locus` is `chr:start-end`
(1-based inclusive) spanning every coordinate of the event — paste it into IGV.
`event_coords` keeps the native 0-based coordinates named per event type:
`exon=…;upstream=…;downstream=…` for SE, `long=…;short=…;flanking=…` for A5SS/A3SS,
`1stExon=…;2ndExon=…;upstream=…;downstream=…` for MXE and
`riExon=…;upstream=…;downstream=…` for RI. The legacy
`exon_start_0base…downstream_ee` columns follow the SE event shape and are empty
for non-SE event types; they are scheduled for removal in 3.0.0.

Per-tool significance rules (identical to the report):

| Tool | Rule |
|------|------|
| rMATS | `padj <= fdr_cutoff & |inc_level_difference| >= dpsi_cutoff` |
| MAJIQ | `probability_changing >= 0.95 & |dpsi_mean| >= 0.2` |
| ISAR | `isoform_switch_q_value < fdr_cutoff & |dIF| >= dpsi_cutoff` |
| LeafCutter | `status == "Success" & padj_cluster <= fdr_cutoff & |deltapsi| >= dpsi_cutoff` (per intron) |
| PEGASAS | no boolean call: exploratory; `n_sig_events` is the BH-adjusted event x pathway correlation count (analytic Pearson p; `n_sig_events_perm` is the permutation cross-check), rank pathways by it |

### cross_tool tables

- `cross_tool.master.tsv` — long format, one row per gene x tool.
- `cross_tool.gene_summary.tsv` — one row per gene: `n_tools_significant`, `tools`,
  `best_effect_size`, `best_padj`. This is the quick gene list to open first.

### QA validator

`run_info/qa_report.txt` is produced by `validate_results_contract.py` after all
deliverables are published, using `assets/results_schema.yaml` as the source of truth for
expected paths and columns. The run aborts if structural pieces are missing (e.g. a master
table for an enabled tool, or an indexed sashimi PDF that does not exist), and also if any
rMATS row has an empty `event_locus` (contract 2.1.0 — events must be localizable in IGV).
Content warnings (FDR==0 fraction, pinned ISAR q-values, tools with 0 significant genes in the cross-tool
table, ...) are reported without blocking.

---

## rMATS Output

**Location**: `results/raw/rmats/<comparison_id>/` (audit-only)

rMATS-turbo detects five types of alternative splicing events and performs statistical testing for differential splicing between conditions.

### Alternative Splicing Event Types

| Event Type | File | Description |
|------------|------|-------------|
| **SE** (Skipped Exon) | `SE.MATS.JC.txt` | Exon inclusion/exclusion (most common AS type) |
| **A5SS** (Alternative 5' Splice Site) | `A5SS.MATS.JC.txt` | Alternative donor site usage |
| **A3SS** (Alternative 3' Splice Site) | `A3SS.MATS.JC.txt` | Alternative acceptor site usage |
| **MXE** (Mutually Exclusive Exons) | `MXE.MATS.JC.txt` | Only one of two exons included |
| **RI** (Retained Intron) | `RI.MATS.JC.txt` | Intron retention events |

### Primary Output Files (per event type)

#### `[AS].MATS.JC.txt` - Junction-Reads-Only Results

**This is the primary output** - uses only reads spanning splice junctions (most reliable).

**Key columns**:

| Column | Description | Example |
|--------|-------------|---------|
| `ID` | Event identifier | `1` |
| `GeneID` | Ensembl gene ID | `ENSG00000123456` |
| `geneSymbol` | Gene name | `BRCA1` |
| `chr` | Chromosome | `chr17` |
| `strand` | + or - | `+` |
| `exonStart_0base` | Event start position | `41234567` |
| `exonEnd` | Event end position | `41234789` |
| `IncLevel1` | Mean PSI in group1 (comma-separated replicates) | `0.85,0.88,0.82` |
| `IncLevel2` | Mean PSI in group2 (comma-separated replicates) | `0.45,0.42,0.48` |
| `IncLevelDifference` | deltaPSI (group2 - group1) | `-0.40` |
| `PValue` | Statistical significance | `0.00012` |
| `FDR` | False Discovery Rate (adjusted p-value) | `0.0034` |

**Interpretation**:
- **PSI** (Percent Spliced In): 0.0 = exon fully excluded, 1.0 = exon fully included
- **deltaPSI > 0**: Increased inclusion in group2 vs group1
- **deltaPSI < 0**: Decreased inclusion in group2 vs group1

**Recommended filtering criteria**:
```r
# Significant events
FDR <= 0.05 & abs(IncLevelDifference) >= 0.1
```

#### `[AS].MATS.JCEC.txt` - Junction + Exon Body Reads

Uses both junction-spanning reads AND exon body reads. Generally more sensitive but can have more false positives. Use for validation/confirmation.

### Annotation Files

| File | Description | Use case |
|------|-------------|----------|
| `fromGTF.[AS].txt` | All AS events detected from GTF + RNA-seq data | Event inventory |
| `fromGTF.novelJunction.[AS].txt` | Novel combinations of known splice sites | Novel transcript isoforms |
| `fromGTF.novelSpliceSite.[AS].txt` | Events with unannotated splice sites | Discovery of new exons (requires `--rmats_novelss true`) |

### Count Files

| File | Description |
|------|-------------|
| `JC.raw.input.[AS].txt` | Raw junction read counts per event |
| `JCEC.raw.input.[AS].txt` | Raw junction + exon counts per event |
| `summary.txt` | Overall statistics: total events, significant events per type |

### Example Usage in R

```r
library(tidyverse)

# Read rMATS results
se_events <- read_tsv("results/raw/rmats/control_vs_treatment/SE.MATS.JC.txt")

# Filter significant events
sig_events <- se_events %>%
  filter(FDR <= 0.05, abs(IncLevelDifference) >= 0.1) %>%
  arrange(FDR)

# Top 10 most differentially spliced events
top10 <- sig_events %>%
  slice_head(n = 10) %>%
  select(geneSymbol, IncLevelDifference, FDR)

print(top10)
```

---

## MAJIQ Output

**Location**: `results/raw/majiq/<comparison_id>/` (audit-only)

MAJIQ uses a Bayesian framework to quantify Local Splicing Variations (LSVs) and compute deltaPSI posteriors.

### Directory Structure

```
majiq/<comparison_id>/
├── splicegraph/
│   ├── <sample1>.majiq
│   ├── <sample2>.majiq
│   └── ...
├── deltapsi/
│   └── <comparison>.deltapsi.voila
└── deltapsi.tsv
```

### Primary Output: `deltapsi.tsv`

**This is the main human-readable output** - TSV export of all LSV results.

**Key columns**:

| Column | Description | Example |
|--------|-------------|---------|
| `gene_name` | Gene symbol | `TP53` |
| `gene_id` | Gene identifier | `ENSG00000141510` |
| `lsv_id` | Local Splicing Variation ID | `s:chr17:7571720:7573927:+` |
| `lsv_type` | Type of variation | `binary` or `complex` |
| `num_junctions` | Number of junctions in LSV | `3` |
| `junctions_coords` | Genomic coordinates | `chr17:7571720-7572927,chr17:7571720-7573927` |
| `de_novo_junctions` | Novel junctions | `0` or `1` (binary) |
| `mean_dpsi_per_lsv_junction` | Expected deltaPSI | `0.35,-0.35` |
| `probability_changing` | P(|\|deltaPSI\|| > threshold) | `0.98` |
| `probability_non_changing` | P(|\|deltaPSI\|| < threshold) | `0.02` |
| `group1_mean_psi` | Mean PSI in group1 | `0.75,0.25` |
| `group2_mean_psi` | Mean PSI in group2 | `0.40,0.60` |

**Interpretation**:
- **LSV** (Local Splicing Variation): A genomic location with alternative junction usage
- **probability_changing ≥ 0.95**: High-confidence differential splicing (default threshold)
- **mean_dpsi_per_lsv_junction**: Comma-separated deltaPSI for each junction in the LSV

**Recommended filtering**:
```python
# Python example
import pandas as pd

majiq = pd.read_csv("results/raw/majiq/control_vs_treatment/deltapsi.tsv", sep="\t")

# High-confidence changing LSVs
sig_lsvs = majiq[majiq['probability_changing'] >= 0.95]

# Further filter by deltaPSI magnitude
high_dpsi = sig_lsvs[
    sig_lsvs['mean_dpsi_per_lsv_junction'].str.extract(r'([\d.]+)').astype(float).max(axis=1) >= 0.2
]
```

### Binary Files (for VOILA visualization)

| File | Description | Tool |
|------|-------------|------|
| `splicegraph/*.majiq` | Per-sample PSI posteriors | Load in MAJIQ VOILA |
| `deltapsi/*.deltapsi.voila` | deltaPSI posteriors | Load in MAJIQ VOILA |
| `splicegraph.sql` | SQLite splice graph database | MAJIQ VOILA or any SQLite browser |

**Visualization**:
```bash
# Use MAJIQ VOILA (separate tool, requires license) to generate interactive visualizations
voila view splicegraph.sql deltapsi/*.deltapsi.voila -o voila_output/
```

---

## IsoformSwitchAnalyzeR Output

**Location**: `results/raw/isar/<comparison_id>/` (audit-only)

IsoformSwitchAnalyzeR identifies isoform switches with predicted functional consequences.

### Primary Outputs

#### `top_isoform_switches.csv`

**Main result table** with all significant isoform switches and their consequences.
The shippable version of this table is `tables/<comparison_id>.isar.master.tsv`
(standard columns + these native ones).

**Key columns**:

| Column | Description | Example |
|--------|-------------|---------|
| `gene_name` | Gene symbol | `MBNL1` |
| `gene_id` | Gene ID | `ENSG00000152601` |
| `isoform_id` | Transcript ID | `ENST00000399503` |
| `condition_1` | Expression in condition 1 | `control` |
| `condition_2` | Expression in condition 2 | `treatment` |
| `IF1` | Isoform Fraction in condition 1 | `0.75` |
| `IF2` | Isoform Fraction in condition 2 | `0.25` |
| `dIF` | Difference in IF (IF2 - IF1) | `-0.50` |
| `isoform_switch_q_value` | Adjusted p-value | `0.0012` |
| `gene_switch_q_value` | Gene-level adjusted p-value | `0.0045` |
| **Consequence columns**: | | |
| `domains_identified` | Protein domain annotation | `PF00018:SH3` |
| `domains_affected` | Domain gain/loss | `gained,lost,none` |
| `signal_peptide_identified` | Signal peptide present | `yes/no` |
| `signal_peptide_affected` | SP gain/loss | `gained/lost/none` |
| `coding_potential` | Coding or noncoding | `coding/noncoding` |
| `ORF_length` | Open reading frame length | `1245` |
| `ORF_seq_similarity` | Sequence similarity between isoforms | `0.87` |
| `NMD_status` | Nonsense-mediated decay sensitivity | `sensitive/insensitive` |

**Interpretation**:
- **IF** (Isoform Fraction): Proportion of gene expression from this isoform (0-1)
- **dIF > 0**: Isoform more abundant in condition 2
- **dIF < 0**: Isoform less abundant in condition 2
- **Switch**: When one isoform increases while another decreases (reciprocal change)

**Filtering example**:
```r
library(tidyverse)

switches <- read_tsv("results/deliverables/contrasts/control_vs_treatment/tables/control_vs_treatment.isar.master.tsv")

# Significant switches with functional consequences
sig_with_consequences <- switches %>%
  filter(isoform_switch_q_value <= 0.05, abs(dIF) >= 0.1) %>%
  filter(domains_affected != "none" | signal_peptide_affected != "none")

# Switches affecting protein domains
domain_switches <- sig_with_consequences %>%
  filter(str_detect(domains_affected, "gained|lost"))
```

#### `consequence_summary.csv`

Summary table counting functional consequence types across all switches
(produced as `consequence_summary.csv` under `raw/isar/<comparison_id>/<comparison_id>/`).

**Columns**:
- `consequence_type`: e.g., "Domain gain", "Signal peptide loss", "NMD sensitive"
- `count`: Number of switches with this consequence
- `proportion`: Fraction of total switches

**Example**:
```
consequence_type          count  proportion
Domain changes              45       0.23
Signal peptide changes      12       0.06
ORF length changes         123       0.63
NMD status changes          34       0.17
```

#### `<comparison_id>_final.rds`

Complete R object containing all data and analysis results. Load in R for custom downstream analyses:

```r
library(IsoformSwitchAnalyzeR)

# Load the switchAnalyzeRlist
aSwitchList <- readRDS("results/raw/isar/control_vs_treatment/control_vs_treatment/control_vs_treatment_final.rds")

# Custom plots
switchPlotTopSwitches(aSwitchList, n = 10, pathToOutput = "my_plots/")

# Extract specific data
isoform_features <- extractSwitchSummary(aSwitchList)
```

#### `switchplots/`

Per-switch PDF plots written by `switchPlotTopN()` (top switches with
consequences), under `raw/isar/<comparison_id>/<comparison_id>/switchplots/`.

---

## LeafCutter Output

**Location**: `results/raw/leafcutter/<comparison_id>/` (audit-only)

LeafCutter quantifies intron usage ratios (not exon inclusion) and is especially sensitive to complex splicing and unannotated introns.

### Files

| File | Description |
|------|-------------|
| `<comparison_id>_cluster_significance.txt` | Per-cluster differential splicing p-values and FDR |
| `<comparison_id>_effect_sizes.txt` | Per-intron log effect size and deltaPSI |

### `<comparison_id>_cluster_significance.txt` columns

| Column | Description |
|--------|-------------|
| `cluster` | Cluster ID (chr:start:end) |
| `p` | Nominal p-value (Chi-squared likelihood ratio test) |
| `p.adjust` | BH-adjusted p-value |
| `df` | Degrees of freedom (introns in cluster − 1) |
| `verdict` | `Significant` or `NotSignificant` |
| `genes` | Associated gene names (when GTF annotation provided) |

### `<comparison_id>_effect_sizes.txt` columns

| Column | Description |
|--------|-------------|
| `intron` | Intron coordinates (chr:start:end:clu_N_strand) |
| `logef` | Log effect size |
| `deltapsi` | Change in intron usage (group2 − group1) |
| `psi1` | Mean intron usage in group1 |
| `psi2` | Mean intron usage in group2 |

---

## Consolidated Report

**Location**: `results/deliverables/contrasts/<comparison_id>/<comparison_id>_splicing_report.html`

### Overview

Interactive HTML report integrating all three tools' results with visualizations and interactive tables.

**Report sections**:

1. **Summary**
   - Experiment overview
   - Sample counts and comparison
   - Tool versions and parameters

2. **rMATS Results**
   - Event counts per AS type (bar chart)
   - Volcano plots (deltaPSI vs -log10(FDR))
   - Interactive data table with top events (with prioritization score: −log10(FDR) × |ΔΨ|)
   - De novo splice-site event table
   - Distribution of deltaPSI values
   - PSI PCA across samples

3. **MAJIQ Results**
   - LSV detection summary
   - deltaPSI distributions (histogram)
   - Probability distributions
   - Interactive table of high-confidence LSVs

4. **IsoformSwitchAnalyzeR Results**
   - Isoform switch counts
   - Functional consequence bar charts
   - Top switches with consequences table
   - Gene-level switch summary

5. **LeafCutter Results** (when `--run_leafcutter true`)
   - Cluster significance table
   - Effect size distribution

6. **PEGASAS Results** (when `--run_pegasas true`)
   - KS score plots per pathway
   - High-correlation event tables

7. **Cross-Tool Overlap**
   - UpSet plot of genes significant across rMATS, MAJIQ, ISAR, and LeafCutter
   - List of high-confidence genes (found by multiple tools)

8. **DE + AS Integration** (when `--de_results` provided)
   - Dual-hit volcano plot overlaying DESeq2/edgeR results with AS hits

9. **GO/KEGG Enrichment**
    - clusterProfiler enrichment for differentially spliced gene sets

10. **Methods**
    - Auto-generated methods section (copy-paste ready for papers)
    - Tool citations
    - Parameter settings

12. **Session Info**
    - R/Python package versions
    - Complete reproducibility information

### Interactive Features

- **DataTables**: Sortable, searchable, exportable result tables
- **Plotly**: Interactive volcano plots with hover labels
- **Collapsible sections**: Show/hide detailed results
- **Download buttons**: Export filtered data as CSV

The report is self-contained: sashimi plots are embedded as PNGs
(resolution controlled by `--sashimi_png_dpi`, default 150) so the HTML opens
offline with no companion files. The vector PDFs are published inside the
deliverables folder under
`deliverables/contrasts/<comparison_id>/plots/sashimi/<EVENT_TYPE>/` and indexed
by `deliverables/contrasts/<comparison_id>/plots/sashimi/sashimi_index.tsv`
(relative `pdf_path`; the `site_class` column is `annotated` or `de_novo` per
event, taken from the SASHIMI_PLOTS manifest `site_classes.tsv` — de novo
events are kept and plotted because `rmats2sashimiplot` can draw them from
rMATS 4.3 coordinates). If the report
container has no `pdftoppm`, sashimi PDFs fall back to being embedded directly
(functional, but heavier).

The report splits the rMATS section into three tabs — Combined (default),
Annotated splice sites and De novo splice sites — each with its own summary,
volcano, top-events ranking, PSI PCA, junction-coverage plot and sashimi
browser. Splice-site class comes from `is_novel_splice_site`
(`fromGTF.novelSpliceSite`): it reflects **splice-site novelty only**, so
events with known splice sites but a novel junction (`fromGTF.novelJunction`)
remain in the annotated tab. Ranking tables filter to events with >=
`report_min_reads` junction reads (IJC+SJC summed per group); significance
calling is never coverage-filtered, and the summary's `Significant_min_reads`
column counts how many significant events clear the minimum. Saturated
priority scores (past the -log10(FDR) cap of 50) are tie-broken by the
weakest group's coverage, and the PSI PCA applies the same coverage filter to
its event pool. The junction-coverage barplot uses a fixed-seed random sample
of SE events with counts (not the first rows in file order).

---

## Upstream QC (MultiQC)

MultiQC is **not** run by this pipeline: it is delivered by the upstream
nf-core/rnaseq run that produces the BAMs. Point `--nfcore_multiqc_dir` (or reuse
that run's `multiqc/multiqc_report.html`) if you need the QC report alongside the
splicing results.

---

## Pipeline Info

**Location**: `results/deliverables/run_info/pipeline_info/` (audit-only)

Nextflow generates these execution artefacts (enabled in `nextflow.config`):

| File | Description |
|------|-------------|
| `execution_report.html` | Visual summary of pipeline execution |
| `execution_timeline.html` | Timeline of process execution (useful for optimization) |
| `execution_trace.txt` | Detailed resource usage per process (CPU, memory, time) |
| `pipeline_dag.svg` | Workflow DAG |

**Use cases**:
- **Debugging**: Identify failed processes
- **Optimization**: Find resource bottlenecks
- **Reporting**: Document compute resources used

---

## File Formats Reference

| Extension | Format | Description | Open with |
|-----------|--------|-------------|-----------|
| `.txt` / `.tsv` | Tab-separated values | Human-readable tables | Excel, R, Python |
| `.csv` | Comma-separated values | Human-readable tables | Excel, R, Python |
| `.rds` | R binary | Serialized R objects | R (`readRDS()`) |
| `.sql` | SQLite database | Relational database | SQLite, DB Browser, MAJIQ VOILA |
| `.majiq` | Binary | MAJIQ PSI posteriors | MAJIQ VOILA |
| `.voila` | Binary | MAJIQ deltaPSI posteriors | MAJIQ VOILA |
| `.html` | HTML | Interactive reports | Any web browser |
| `.fasta` | FASTA | Sequence data | Text editor, bioinformatics tools |

---

## Data Retention and Cleanup

### Intermediate Files (`work/` directory)

The `work/` directory contains all intermediate files from Nextflow processes:
- Temporary BAM subsets
- rMATS prep files
- Intermediate R objects

**Size**: Can be 10-50× larger than `results/`

**Retention**: 
- Keep if you plan to re-run with `-resume`
- Delete after successful completion to free space:
  ```bash
  rm -rf work/
  # Or use Nextflow's cleanup
  nextflow clean -f -k
  ```

### Long-term Storage

Only `results/deliverables/` ships; `results/raw/` is audit evidence. What to keep:

- **Ship (compact):** the whole `deliverables/` folder (rename it and zip it).
  Optional: drop `run_info/pipeline_info/` if you do not need the Nextflow trace.
- **Minimal set** (for publications): each `deliverables/contrasts/<comparison>/`
  folder — the HTML report, `tables/*.master.tsv` / `*.significant.tsv` and
  `plots/sashimi/`.
- **Audit (keep, do not ship):** `results/raw/` (native tool outputs, MAJIQ
  binaries, rMATS fromGTF/JCEC, sashimi intermediates), `results/run_info/`
  (`run_manifest.yaml` + per-contrast manifests) and
  `deliverables/run_info/pipeline_info/`.

---

## Example Analysis Workflows

### Read the standardised tables (recommended)

```r
library(tidyverse)

comp <- "control_vs_treatment"
tables <- sprintf("results/deliverables/contrasts/%s/tables", comp)

# One row per gene with the tools that called it significant
gene_summary <- read_tsv(file.path(tables, sprintf("%s.cross_tool.gene_summary.tsv", comp)))
top_genes <- gene_summary %>% filter(n_tools_significant >= 2)

# Full rMATS master (all events), then the significant subset
rmats <- read_tsv(file.path(tables, sprintf("%s.rmats.master.tsv", comp)))
rmats_sig <- read_tsv(file.path(tables, sprintf("%s.rmats.significant.tsv", comp)))
```

`significant.tsv` already encodes the report rule; use `padj` + `padj_method`
rather than a raw `fdr` column.

### Export a gene list for enrichment

```python
import pandas as pd

tables = "results/deliverables/contrasts/control_vs_treatment/tables"
sig = pd.read_csv(f"{tables}/control_vs_treatment.rmats.significant.tsv", sep="\t")
gene_list = sig["gene_symbol"].dropna().unique().tolist()
print(len(gene_list), "genes")
```

---

## Citation

When publishing results, cite all tools used:

- **rMATS**: Shen et al. (2014) PNAS. doi:10.1073/pnas.1419161111
- **MAJIQ**: Vaquero-Garcia et al. (2016) eLife. doi:10.7554/eLife.11752
- **IsoformSwitchAnalyzeR**: Vitting-Seerup & Sandelin (2019) Mol Cell. doi:10.1016/j.molcel.2019.09.005

Full citations available in `CITATIONS.md`.
