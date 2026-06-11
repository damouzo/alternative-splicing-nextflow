# Usage Documentation

## Introduction

This pipeline analyzes differential alternative splicing from pre-processed RNA-seq data produced by [nf-core/rnaseq](https://nf-co.re/rnaseq).

## Prerequisites

### Upstream nf-core/rnaseq Run

This pipeline requires outputs from nf-core/rnaseq v3.14.0 or later, run with:
- `--aligner star_salmon` (default)
- `--save_align_intermeds` (to preserve BAM files)
- `--save_quantification_intermeds` (for Salmon directories)

Optionally, Salmon should be run with `--numBootstraps 100` for improved IsoformSwitchAnalyzeR accuracy.

### Required Files Per Sample

From the nf-core/rnaseq `star_salmon/<sample>/` directory:
- Coordinate-sorted BAM: `<sample>.markdup.sorted.bam`
- BAM index: `<sample>.markdup.sorted.bam.bai`
- Salmon quantification directory containing `quant.sf`

### Reference Files

1. **GTF annotation** - MUST be the exact same file used in the nf-core/rnaseq run
2. **Genome FASTA** - Required when using `--use_gffread true` (recommended)
3. **Transcript FASTA** - Optional when using gffread; otherwise must match GTF release exactly

## Running the Pipeline

### Step 1: Prepare Samplesheet

Create a CSV file with the following columns:

```csv
sample,condition,replicate,bam,bai,salmon_dir
ctrl_rep1,control,1,/path/to/ctrl_rep1.markdup.sorted.bam,/path/to/ctrl_rep1.markdup.sorted.bam.bai,/path/to/ctrl_rep1/
ctrl_rep2,control,2,/path/to/ctrl_rep2.markdup.sorted.bam,/path/to/ctrl_rep2.markdup.sorted.bam.bai,/path/to/ctrl_rep2/
ctrl_rep3,control,3,/path/to/ctrl_rep3.markdup.sorted.bam,/path/to/ctrl_rep3.markdup.sorted.bam.bai,/path/to/ctrl_rep3/
treat_rep1,treatment,1,/path/to/treat_rep1.markdup.sorted.bam,/path/to/treat_rep1.markdup.sorted.bam.bai,/path/to/treat_rep1/
treat_rep2,treatment,2,/path/to/treat_rep2.markdup.sorted.bam,/path/to/treat_rep2.markdup.sorted.bam.bai,/path/to/treat_rep2/
treat_rep3,treatment,3,/path/to/treat_rep3.markdup.sorted.bam,/path/to/treat_rep3.markdup.sorted.bam.bai,/path/to/treat_rep3/
```

**Important requirements:**
- Minimum 3 replicates per condition (recommended for statistical power)
- Exactly 2 conditions
- All file paths must exist and be readable

### Step 2: Prepare Comparisons File

Create a CSV specifying which condition is group1 (reference) and group2 (treatment):

```csv
group1,group2
control,treatment
```

### Step 3: Determine Required Parameters

Check your nf-core/rnaseq MultiQC report for:

1. **Strandedness**: Look in the "Strand specificity" section or check any sample's `lib_format_counts.json`
   - `unstranded`, `forward`, or `reverse`

2. **Read Length**: Post-trimming read length from alignment summary
   - The pipeline will verify this against actual BAM content

### Step 4: Run the Pipeline

Using a parameters file (recommended):

```bash
# Create params.yaml
cat > params.yaml <<'EOF'
input: samplesheet.csv
comparisons: comparisons.csv
outdir: ./results

# References (must match nf-core/rnaseq references)
gtf: /ref/Homo_sapiens.GRCh38.110.gtf
genome_fasta: /ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa
use_gffread: true
genome_build: hg38

# Library properties
strandedness: reverse
read_length: 150

# Tool switches
run_rmats: true
run_majiq: true
run_isar: true

# MAJIQ license (required)
majiq_license: /path/to/majiq_license.key
EOF

# Run pipeline
nextflow run main.nf \
  -profile docker \
  -params-file params.yaml
```

Or using command-line arguments:

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --comparisons comparisons.csv \
  --gtf reference.gtf \
  --genome_fasta genome.fa \
  --use_gffread true \
  --strandedness reverse \
  --read_length 150 \
  --majiq_license majiq_license.key \
  --outdir results \
  -profile docker
```

### Step 5: Monitor Execution

Nextflow will display progress in real-time:

```
executor >  local (12)
[4a/b5c123] process > ALTERNATIVE_SPLICING:INPUT_CHECK:VALIDATE_INPUT       [100%] 1 of 1 ✔
[2f/8e9abc] process > ALTERNATIVE_SPLICING:RMATS_ANALYSIS:RMATS_PREP (1)    [ 50%] 3 of 6
[7d/123def] process > ALTERNATIVE_SPLICING:MAJIQ_ANALYSIS:MAJIQ_BUILD       [  0%] 0 of 1
```

## Output Structure

Results are written to the directory specified by `--outdir`:

```
results/
├── rmats/
│   └── control_vs_treatment/
│       ├── SE.MATS.JC.txt
│       ├── A5SS.MATS.JC.txt
│       └── ...
├── majiq/
│   └── control_vs_treatment/
│       ├── splicegraph.sql
│       ├── control_vs_treatment.tsv
│       └── ...
├── isoformswitchr/
│   └── control_vs_treatment/
│       ├── top_isoform_switches.csv
│       ├── consequence_summary.csv
│       └── switchplots/
└── report/
    └── control_vs_treatment_splicing_report.html
```

See [output.md](output.md) for detailed descriptions of all output files.

## Common Issues

### Chromosome Naming Mismatch

**Error**: "No events detected" or "0 splice junctions found" or empty MAJIQ splice graph

**Cause**: GTF uses `chr1, chr2, ...` but BAM uses `1, 2, ...` (or vice versa)

**Detection**: Pipeline validates this during `INPUT_CHECK` and will fail early with:
```
ERROR: Chromosome naming mismatch detected
  BAM chromosomes: 1, 2, 3, ...
  GTF chromosomes: chr1, chr2, chr3, ...
```

**Solution**: Use the **exact same GTF** in both nf-core/rnaseq and this pipeline. The alignment and annotation must use consistent naming.

---

### Read Length Mismatch

**Error**: "Read length validation failed: expected 150, found 140 in sample X"

**Cause**: `--read_length` doesn't match actual post-trimming read length in BAMs

**Why it matters**: rMATS uses read length for junction detection sensitivity, MAJIQ uses it in `readlen` config parameter. Incorrect values cause biased PSI estimates.

**Solution**: 
1. Check actual read length:
   ```bash
   samtools view sample.bam | head -1000 | awk '{print length($10)}' | sort | uniq -c
   ```
2. Or check nf-core/rnaseq MultiQC report → "Alignment" section
3. Update `--read_length` to match the modal value

---

### Insufficient Replicates

**Warning**: "Condition 'control' has only 2 replicates (minimum 3 recommended)"

**Impact**: 
- Reduced statistical power
- Unreliable FDR estimates
- High false discovery rate
- Cannot distinguish biological variation from technical noise

**Solution**: 
- Ideal: Add more biological replicates (3-5 per condition)
- Workaround: Use more lenient significance thresholds (e.g., FDR < 0.1) but report as exploratory
- **Do not**: Pool technical replicates as biological replicates

---

### MAJIQ License Missing

**Error**: "MAJIQ_LICENSE_FILE not set or invalid"

**Solution**: 
1. Obtain academic license from https://majiq.biociphers.org (free for academic use)
2. Pass via parameter:
   ```bash
   nextflow run main.nf --majiq_license /path/to/license.key ...
   ```
3. Or set environment variable:
   ```bash
   export MAJIQ_LICENSE_FILE=/path/to/license.key
   ```

---

### Out of Memory Errors

**Error**: "Process exceeded available memory"

**Symptoms**:
- Process killed with exit code 137 (OOM killed)
- Java heap space errors
- Segmentation faults in C++ tools

**Solutions**:

1. **Increase memory limits**:
   ```bash
   nextflow run main.nf --max_memory '256.GB' ...
   ```

2. **Reduce parallelism** (process fewer samples simultaneously):
   Edit `nextflow.config`:
   ```groovy
   executor {
       queueSize = 4  // Reduce from default
   }
   ```

3. **For Nextflow JVM** (not tool processes):
   ```bash
   export NXF_OPTS='-Xms2g -Xmx8g'
   ```

4. **Memory-intensive processes**: rMATS POST and MAJIQ BUILD typically need 32-64 GB for human genome

---

### No Significant Events Detected

**Symptom**: All output tables are empty or have 0 significant events

**Possible causes**:

1. **True biological result**: No differential splicing between conditions
   - Check: Are conditions actually different? (drug vs vehicle, tissue A vs B)
   - Verify: Check gene-level differential expression (nf-core/rnaseq results)

2. **Insufficient coverage**: Low read depth at junctions
   - Check: `multiqc_report.html` → "Alignment metrics" → junction counts
   - Solution: Sequence deeper (aim for 50-100M reads/sample)

3. **Poor RNA quality**: 3' bias masking internal splice junctions
   - Check: RSeQC gene body coverage plots (should be uniform)
   - Solution: Re-do library prep with higher quality RNA (RIN ≥ 7)

4. **Wrong library strandedness**: Reduced junction detection
   - Check: Salmon `lib_format_counts.json` → should match `--strandedness`
   - Solution: Correct `--strandedness` parameter

5. **Overly stringent thresholds**:
   - Try relaxing: `--rmats_cutoff 0.1` or `--majiq_deltapsi_threshold 0.1`

---

### Docker/Singularity Issues

**Error**: "Cannot pull image" or "Singularity not found"

**Solutions**:

1. **Docker not running**:
   ```bash
   sudo systemctl start docker
   ```

2. **Insufficient disk space** for images:
   ```bash
   docker system prune -a  # Clean up old images
   ```

3. **Singularity cache full**:
   ```bash
   export NXF_SINGULARITY_CACHEDIR=/path/to/large/disk
   ```

4. **Behind firewall** (cannot pull images):
   - Pre-pull all images manually
   - Or build containers locally from Dockerfiles in `containers/`

---

### Pipeline Hangs or Stalls

**Symptoms**: Process shows 0% progress for hours, no CPU/memory usage

**Causes**:

1. **Waiting for cluster resources**: Check queue status
   ```bash
   squeue -u $USER  # Slurm
   qstat -u $USER   # SGE
   ```

2. **Input file locked**: Check if BAM/BAI files are being written by another process

3. **Network filesystem latency**: Temporary NFS hang
   - Wait or restart Nextflow with `-resume`

4. **Java garbage collection**: Nextflow JVM paused
   - Increase heap size: `export NXF_OPTS='-Xmx8g'`

---

### Comparison ID Not Found

**Error**: "Sample X has condition 'ctrl' which is not in comparisons.csv"

**Cause**: Mismatch between `samplesheet.csv` conditions and `comparisons.csv` groups

**Solution**: Ensure every `condition` value in samplesheet appears in either `group1` or `group2` column of comparisons.csv

**Example**:
```csv
# samplesheet.csv
sample,condition,replicate,...
s1,control,1,...
s2,treatment,1,...

# comparisons.csv
group1,group2
control,treatment  # ✓ Both conditions present
```

---

### IsoformSwitchAnalyzeR Errors

**Error**: "No isoforms pass filtering criteria"

**Causes**:
1. **Low Salmon bootstrap count**: Re-run nf-core/rnaseq with `--salmon_bootstraps 100`
2. **Insufficient isoform expression**: Most transcripts have IF < 0.01
3. **No transcript FASTA**: Check `--transcript_fasta` or `--use_gffread true` is set

**Error**: "GTF and transcript FASTA mismatch"

**Cause**: FASTA from different Ensembl/Gencode release than GTF

**Solution**: Download matching release:
```bash
# Gencode v44 example
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
```

Or use `--use_gffread true` to extract from genome+GTF (ensures consistency).

## Complete Parameter Reference

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--input` | file | Path to samplesheet CSV (see format above) |
| `--comparisons` | file | Path to comparisons CSV (group1,group2) |
| `--gtf` | file | Gene annotation GTF (must match nf-core/rnaseq reference) |
| `--strandedness` | string | Library strandedness: `unstranded`, `forward`, or `reverse` |
| `--read_length` | integer | Post-trimming read length (validated against BAMs) |

### Reference Files

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--genome_fasta` | file | `null` | Genome FASTA (required if `use_gffread=true`) |
| `--transcript_fasta` | file | `null` | Transcript FASTA (alternative to gffread, must match GTF) |
| `--use_gffread` | boolean | `true` | Extract transcriptome from genome+GTF using gffread |
| `--genome_build` | string | `'hg38'` | Genome build identifier for MAJIQ (e.g., 'hg38', 'mm10', 'GRCh38') |

**Note**: Either provide `--genome_fasta` (with `use_gffread=true`) OR provide `--transcript_fasta` (with `use_gffread=false`). The transcript FASTA must exactly match the GTF gene model release.

### Tool Toggles

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--run_rmats` | boolean | `true` | Run rMATS-turbo differential splicing analysis |
| `--run_majiq` | boolean | `true` | Run MAJIQ deltaPSI quantification |
| `--run_isar` | boolean | `true` | Run IsoformSwitchAnalyzeR isoform switch analysis |
| `--run_sashimi` | boolean | `false` | Generate sashimi plots for top rMATS events |
| `--run_pegasas` | boolean | `false` | Run PEGASAS pathway–splicing PSI correlation |
| `--run_leafcutter` | boolean | `false` | Run LeafCutter intron excision analysis |
| `--run_isar_full_annotation` | boolean | `false` | Enable Tier A ISAR annotation (PFAM + IUPred3) |

**Example**: To run only rMATS and skip MAJIQ/ISAR:
```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --run_rmats true \
  --run_majiq false \
  --run_isar false \
  ...
```

### rMATS Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--rmats_novel_ss` | boolean | `false` | Enable novel splice site detection (increases runtime ~3-5×) |
| `--rmats_cstat` | float | `0.0001` | rMATS Cstat significance threshold |
| `--rmats_min_intron_length` | integer | `50` | Minimum intron length |
| `--rmats_max_exon_length` | integer | `500` | Maximum exon length for filtering |
| `--rmats_tstat_threads` | integer | `6` | Threads for t-statistic computation |
| `--rmats_read_type` | string | `'paired'` | Read type: `paired` or `single` |

### MAJIQ Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--majiq_license` | file | env `MAJIQ_LICENSE_FILE` | MAJIQ academic licence key file |
| `--majiq_sif` | string | env `MAJIQ_SIF` | MAJIQ container image path/tag |
| `--majiq_min_denovo_reads` | integer | `2` | Min reads supporting a de novo junction |
| `--majiq_min_intronic_cov` | integer | `3` | Min intronic coverage |
| `--majiq_min_reads` | integer | `10` | Min reads for junction detection |
| `--majiq_min_nonzero` | integer | `3` | Min non-zero samples |
| `--majiq_delta_psi_threshold` | float | `0.2` | deltaPSI reporting threshold |
| `--majiq_probability_threshold` | float | `0.95` | Posterior probability threshold |

**Note**: MAJIQ requires an academic license. Set `MAJIQ_LICENSE_FILE` environment variable or pass `--majiq_license`.

### IsoformSwitchAnalyzeR Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--isar_alpha` | float | `0.05` | Significance threshold for isoform switches |
| `--isar_dif_cutoff` | float | `0.1` | Minimum isoform fraction difference |
| `--isar_gene_expr_cutoff` | float | `1` | Minimum gene TPM |
| `--isar_iso_expr_cutoff` | float | `1` | Minimum isoform TPM |
| `--run_isar_full_annotation` | boolean | `false` | Enable Tier A annotation: PFAM + IUPred3 IDR |
| `--pfam_hmm` | file | `null` | Path to Pfam-A.hmm database file (~500 MB) |

### Sashimi Plots

Requires `--run_sashimi true`. Generates arc plots via rmats2sashimiplot for top rMATS SE events.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--run_sashimi` | boolean | `false` | Enable sashimi plot generation |
| `--sashimi_top_n` | integer | `10` | Number of top events to plot |
| `--sashimi_group1_label` | string | `'Group1'` | Label for group1 in plots |
| `--sashimi_group2_label` | string | `'Group2'` | Label for group2 in plots |
| `--sashimi_exon_scale` | integer | `25` | Exon scale factor |
| `--sashimi_intron_scale` | integer | `5` | Intron scale factor |

```bash
nextflow run main.nf -profile docker \
  --run_sashimi true --sashimi_top_n 20 \
  --sashimi_group1_label control --sashimi_group2_label treatment \
  -params-file params.yaml
```

### PEGASAS

Requires `--run_pegasas true`. Correlates pathway activity scores (from gene sets) with per-sample PSI values.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--run_pegasas` | boolean | `false` | Enable PEGASAS pathway–splicing correlation |
| `--pathway_gmt` | file | `null` | GMT file with gene sets (e.g., MSigDB Hallmarks) |
| `--salmon_merged_tpm` | file | `null` | Merged TPM matrix from nf-core/rnaseq (all samples) |
| `--pegasas_groups` | string | `null` | Comma-separated condition labels matching samplesheet |
| `--pegasas_num_interval` | integer | `100` | Number of KS enrichment intervals |

```bash
nextflow run main.nf -profile docker \
  --run_pegasas true \
  --pathway_gmt /data/h.all.v2023.1.Hs.symbols.gmt \
  --salmon_merged_tpm /data/salmon.merged.gene_tpm.tsv \
  -params-file params.yaml
```

### DE + AS Integration

Overlay DESeq2/edgeR results on the splicing report as a dual-hit volcano plot.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--de_results` | file | `null` | TSV with DE results (columns: gene_name, log2FoldChange, padj) |

### LeafCutter

Requires `--run_leafcutter true`. Performs intron excision clustering and differential splicing via LeafCutter.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--run_leafcutter` | boolean | `false` | Enable LeafCutter analysis |
| `--leafcutter_container` | string | GHCR tag | Override LeafCutter container image |

```bash
nextflow run main.nf -profile docker \
  --run_leafcutter true \
  -params-file params.yaml
```

### QC / Report

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--nfcore_multiqc_dir` | path | `null` | Path to nf-core/rnaseq MultiQC output directory |
| `--organism` | string | `'human'` | Organism for GO/KEGG enrichment (`human` or `mouse`) |

### Container Overrides

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--isar_container` | string | GHCR tag | Override ISAR container image |
| `--report_container` | string | GHCR tag | Override report container image |
| `--pegasas_container` | string | GHCR tag | Override PEGASAS container image |
| `--leafcutter_container` | string | GHCR tag | Override LeafCutter container image |

### Output Options

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--outdir` | path | `'./results'` | Output directory for all results |
| `--publish_dir_mode` | string | `'copy'` | File publishing mode: `copy`, `symlink`, or `move` |

### Execution Profiles

Use `-profile <name>` to select execution environment:

| Profile | Description |
|---------|-------------|
| `docker` | Run using Docker containers (recommended for local) |
| `singularity` | Run using Singularity/Apptainer containers (HPC) |
| `apptainer` | Run using Apptainer containers (modern Singularity successor) |
| `slurm` | Submit jobs to a SLURM cluster |
| `lowmem` | Cap ISAR process concurrency to 1 and limit executor resources to 4 CPUs / 16 GB. Use on small workstations where the default uncapped ISAR parallelisation causes OOM. |
| `test` | Run with minimal test data (for pipeline validation) |
| `test_full` | Run with the full demo dataset |

Profiles can be combined. **Example**:
```bash
nextflow run main.nf -profile singularity,slurm -params-file params.yaml

# Low-memory workstation
nextflow run main.nf -profile docker,lowmem -params-file params.yaml
```

## Best Practices

### Nextflow JVM Memory

Prevent out-of-memory errors in Nextflow's Java runtime by setting:

```bash
# Add to ~/.bashrc or export before running
export NXF_OPTS='-Xms1g -Xmx4g'
```

### Reproducibility

Always specify pipeline version when running on production data:

```bash
# Pin to specific release
nextflow run alternative-splicing-nextflow -r v1.0.0 -params-file params.yaml

# Or use a specific git commit
nextflow run alternative-splicing-nextflow -r a1b2c3d -params-file params.yaml
```

### Parameter Files vs Command Line

For complex runs, use YAML parameter files:

**Advantages**:
- Version-controlled alongside code
- Easier to review and share
- Prevents command-line typos
- Documents exact parameters used

```yaml
# params.yaml
input: /data/project/samplesheet.csv
comparisons: /data/project/comparisons.csv
outdir: /data/project/results_v1

gtf: /ref/gencode.v44.annotation.gtf
genome_fasta: /ref/GRCh38.primary_assembly.fa
use_gffread: true
genome_build: 'GRCh38'

strandedness: reverse
read_length: 150

run_rmats: true
run_majiq: true
run_isar: true

majiq_license: /home/user/.majiq/license.key

max_cpus: 32
max_memory: '256.GB'
```

Run with:
```bash
nextflow run main.nf -r v1.0.0 -profile singularity -params-file params.yaml
```

### Recommended Minimum Requirements

For robust differential splicing analysis:

| Requirement | Minimum | Recommended | Rationale |
|-------------|---------|-------------|-----------|
| **Replicates per group** | 2 | 3-5 | Statistical power for FDR control |
| **Read depth** | 20M | 50-100M | Junction coverage for rare events |
| **Read length** | 75 bp | 100-150 bp | Splice junction mappability |
| **CPUs** | 4 | 16-32 | Parallel processing of samples |
| **Memory** | 32 GB | 64-128 GB | rMATS and MAJIQ memory usage |
| **Storage** | 50 GB | 200+ GB | Intermediate files and results |

### Execution Time Estimates

Typical runtime for 6 samples (3 vs 3), human genome, 50M reads/sample:

| Tool | Runtime | Bottleneck |
|------|---------|------------|
| **Input validation** | ~2-5 min | BAM chromosome extraction |
| **rMATS PREP** | ~30-60 min | Per-sample junction extraction |
| **rMATS POST** | ~15-30 min | Statistical testing |
| **MAJIQ BUILD** | ~20-40 min | Splice graph construction |
| **MAJIQ DELTAPSI** | ~10-20 min | Posterior probability computation |
| **ISAR Import** | ~5-10 min | Salmon data loading |
| **ISAR Analysis** | ~20-40 min | ORF extraction and annotation |
| **Report generation** | ~5-10 min | R Markdown rendering |
| **Total** | ~2-4 hours | (with parallelization) |

**Note**: Runtime scales approximately linearly with:
- Number of samples (rMATS PREP, MAJIQ BUILD)
- Genome size (all tools)
- Junction complexity (rMATS POST, MAJIQ DELTAPSI)

## Advanced Options

### Running Subset of Tools

```bash
# Only rMATS (fastest)
nextflow run main.nf --run_rmats true --run_majiq false --run_isar false -params-file params.yaml

# rMATS + MAJIQ (skip ISAR)
nextflow run main.nf --run_isar false -params-file params.yaml
```

### Custom Resource Allocation

Override default resources per-process in `nextflow.config` or via CLI:

```bash
nextflow run main.nf \
  --max_cpus 64 \
  --max_memory '512.GB' \
  -params-file params.yaml
```

### Resume Failed Runs

Nextflow caches completed tasks. If a run fails, fix the issue and resume:

```bash
nextflow run main.nf -resume -params-file params.yaml
```

Only failed/pending tasks will re-run.

### Work Directory Cleanup

After successful completion, clean up intermediate files:

```bash
# Remove work directory (keeps only published results)
rm -rf work/

# Or use Nextflow's built-in cleanup
nextflow clean -f
```

**Warning**: This deletes all intermediate files. Don't do this if you might need to `-resume`.

## Troubleshooting

For complete parameter documentation, run:

```bash
nextflow run main.nf --help
```

Or see the full parameter list in `nextflow.config`.
