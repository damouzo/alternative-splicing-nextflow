# Output Files Documentation

This document describes all output files produced by the alternative-splicing-nextflow pipeline.

## Directory Structure

```
results/
├── rmats/
│   └── <comparison_id>/
│       ├── SE.MATS.JC.txt
│       ├── A5SS.MATS.JC.txt
│       ├── A3SS.MATS.JC.txt
│       ├── MXE.MATS.JC.txt
│       ├── RI.MATS.JC.txt
│       └── summary.txt
├── majiq/
│   └── <comparison_id>/
│       ├── splicegraph/
│       │   └── *.majiq
│       ├── deltapsi/
│       │   └── *.deltapsi.voila
│       └── deltapsi.tsv
├── isoformswitchanalyzer/
│   └── <comparison_id>/
│       ├── switchAnalyzeRlist_final.rds
│       ├── isoform_switches.txt
│       ├── consequences_summary.txt
│       └── isoformSwitchAnalyzeR_AA.fasta
├── report/
│   └── <comparison_id>_splicing_report.html
├── multiqc/
│   └── multiqc_report.html
└── pipeline_info/
    ├── execution_report.html
    ├── execution_timeline.html
    └── execution_trace.txt
```

---

## rMATS Output

**Location**: `results/rmats/<comparison_id>/`

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
se_events <- read_tsv("results/rmats/control_vs_treatment/SE.MATS.JC.txt")

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

**Location**: `results/majiq/<comparison_id>/`

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

majiq = pd.read_csv("results/majiq/control_vs_treatment/deltapsi.tsv", sep="\t")

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

**Location**: `results/isoformswitchanalyzer/<comparison_id>/`

IsoformSwitchAnalyzeR identifies isoform switches with predicted functional consequences.

### Primary Outputs

#### `isoform_switches.txt`

**Main result table** with all significant isoform switches and their consequences.

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

switches <- read_tsv("results/isoformswitchanalyzer/control_vs_treatment/isoform_switches.txt")

# Significant switches with functional consequences
sig_with_consequences <- switches %>%
  filter(isoform_switch_q_value <= 0.05, abs(dIF) >= 0.1) %>%
  filter(domains_affected != "none" | signal_peptide_affected != "none")

# Switches affecting protein domains
domain_switches <- sig_with_consequences %>%
  filter(str_detect(domains_affected, "gained|lost"))
```

#### `consequences_summary.txt`

Summary table counting functional consequence types across all switches.

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

#### `switchAnalyzeRlist_final.rds`

Complete R object containing all data and analysis results. Load in R for custom downstream analyses:

```r
library(IsoformSwitchAnalyzeR)

# Load the switchAnalyzeRlist
aSwitchList <- readRDS("results/isoformswitchanalyzer/control_vs_treatment/switchAnalyzeRlist_final.rds")

# Custom plots
switchPlotTopSwitches(aSwitchList, n = 10, pathToOutput = "my_plots/")

# Extract specific data
isoform_features <- extractSwitchSummary(aSwitchList)
```

#### `isoformSwitchAnalyzeR_AA.fasta`

Amino acid sequences for all annotated isoforms. Use for:
- External domain prediction tools
- Sequence alignment
- Structural modeling

---

## Consolidated Report

**Location**: `results/report/<comparison_id>_splicing_report.html`

### Overview

Interactive HTML report integrating all three tools' results with visualizations and interactive tables.

**Report sections**:

1. **Summary**
   - Experiment overview
   - Sample counts and comparison
   - Tool versions and parameters

2. **Quality Control**
   - Embedded MultiQC metrics (if available)
   - Sample validation results

3. **rMATS Results**
   - Event counts per AS type (bar chart)
   - Volcano plots (deltaPSI vs -log10(FDR))
   - Interactive data table with top events
   - Distribution of deltaPSI values

4. **MAJIQ Results**
   - LSV detection summary
   - deltaPSI distributions (histogram)
   - Probability distributions
   - Interactive table of high-confidence LSVs

5. **IsoformSwitchAnalyzeR Results**
   - Isoform switch counts
   - Functional consequence bar charts
   - Top switches with consequences table
   - Gene-level switch summary

6. **Cross-Tool Overlap**
   - Venn diagram of genes significant in ≥2 tools
   - List of high-confidence genes (found by all 3 tools)

7. **Methods**
   - Auto-generated methods section (copy-paste ready for papers)
   - Tool citations
   - Parameter settings

8. **Session Info**
   - R/Python package versions
   - Complete reproducibility information

### Interactive Features

- **DataTables**: Sortable, searchable, exportable result tables
- **Plotly**: Interactive volcano plots with hover labels
- **Collapsible sections**: Show/hide detailed results
- **Download buttons**: Export filtered data as CSV

---

## MultiQC Report

**Location**: `results/multiqc/multiqc_report.html`

Aggregates QC metrics from all tools into a single report.

**Includes**:
- rMATS event detection statistics
- Sample-level QC metrics (if upstream logs provided)
- Junction read counts
- Splice graph complexity

---

## Pipeline Info

**Location**: `results/pipeline_info/`

Nextflow automatically generates execution reports:

| File | Description |
|------|-------------|
| `execution_report.html` | Visual summary of pipeline execution |
| `execution_timeline.html` | Timeline of process execution (useful for optimization) |
| `execution_trace.txt` | Detailed resource usage per process (CPU, memory, time) |

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

Only `results/` directory needs long-term storage. Key files for publications/downstream analysis:

**Minimal set** (for publications):
- `report/<comparison>_splicing_report.html`
- `rmats/<comparison>/SE.MATS.JC.txt` (and other event types)
- `majiq/<comparison>/deltapsi.tsv`
- `isoformswitchanalyzer/<comparison>/isoform_switches.txt`

**Complete set** (for re-analysis):
- Everything in `results/` (includes binary files for MAJIQ VOILA, R objects, etc.)

---

## Example Analysis Workflows

### Extract Top Events from All Tools

```r
library(tidyverse)

# rMATS top SE events
rmats_se <- read_tsv("results/rmats/control_vs_treatment/SE.MATS.JC.txt") %>%
  filter(FDR <= 0.05, abs(IncLevelDifference) >= 0.1) %>%
  select(geneSymbol, IncLevelDifference, FDR)

# MAJIQ top LSVs
majiq <- read_tsv("results/majiq/control_vs_treatment/deltapsi.tsv") %>%
  filter(probability_changing >= 0.95)

# ISAR top switches
isar <- read_tsv("results/isoformswitchanalyzer/control_vs_treatment/isoform_switches.txt") %>%
  filter(isoform_switch_q_value <= 0.05, abs(dIF) >= 0.1)

# Find genes in common
common_genes <- reduce(
  list(
    rmats = unique(rmats_se$geneSymbol),
    majiq = unique(majiq$gene_name),
    isar = unique(isar$gene_name)
  ),
  intersect
)

print(paste("Genes significant in all 3 tools:", length(common_genes)))
```

### Export for Pathway Analysis

```python
import pandas as pd

# Load rMATS SE results
se = pd.read_csv("results/rmats/control_vs_treatment/SE.MATS.JC.txt", sep="\t")

# Filter significant
sig = se[(se['FDR'] <= 0.05) & (abs(se['IncLevelDifference']) >= 0.1)]

# Extract gene list for GSEA/GO enrichment
gene_list = sig['geneSymbol'].unique().tolist()

# Save for enrichment analysis
with open("genes_with_differential_splicing.txt", "w") as f:
    f.write("\n".join(gene_list))
```

---

## Citation

When publishing results, cite all tools used:

- **rMATS**: Shen et al. (2014) PNAS. doi:10.1073/pnas.1419161111
- **MAJIQ**: Vaquero-Garcia et al. (2016) eLife. doi:10.7554/eLife.11752
- **IsoformSwitchAnalyzeR**: Vitting-Seerup & Sandelin (2019) Mol Cell. doi:10.1016/j.molcel.2019.09.005

Full citations available in `CITATIONS.md`.
