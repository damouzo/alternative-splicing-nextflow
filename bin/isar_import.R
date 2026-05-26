#!/usr/bin/env Rscript

# Import Salmon quantifications into IsoformSwitchAnalyzeR

suppressPackageStartupMessages({
    library(IsoformSwitchAnalyzeR)
})

# Parse --key value arguments from command line
.parse_args <- function() {
    raw <- commandArgs(trailingOnly = TRUE)
    out <- list()
    i <- 1L
    while (i <= length(raw)) {
        if (startsWith(raw[i], "--")) {
            key <- sub("^--", "", raw[i])
            if (i + 1L <= length(raw) && !startsWith(raw[i + 1L], "--")) {
                out[[key]] <- raw[i + 1L]; i <- i + 2L
            } else {
                out[[key]] <- TRUE; i <- i + 1L
            }
        } else { i <- i + 1L }
    }
    out
}
opt <- .parse_args()

# Validate inputs
if (is.null(opt$samplesheet) || is.null(opt$gtf) ||
    is.null(opt$transcript_fasta) || is.null(opt$output)) {
    stop("Required arguments: --samplesheet --gtf --transcript_fasta --output")
}

cat("===========================================\n")
cat("IsoformSwitchAnalyzeR - Import\n")
cat("===========================================\n")

# Read samplesheet
samplesheet <- read.csv(opt$samplesheet, stringsAsFactors = FALSE)
cat("Samples:", nrow(samplesheet), "\n")
cat("Conditions:", length(unique(samplesheet$condition)), "\n")

# Prepare design matrix
design <- data.frame(
    sampleID  = samplesheet$sample,
    condition = samplesheet$condition,
    stringsAsFactors = FALSE
)

# Prepare Salmon quantification paths
salmon_quant_paths <- data.frame(
    sampleID = samplesheet$sample,
    quant_path = file.path(samplesheet$salmon_dir, "quant.sf"),
    stringsAsFactors = FALSE
)

cat("\nImporting Salmon quantifications...\n")

# ISA v2+ API: step 1 — import expression only (no annotation args here)
quant_vec <- setNames(salmon_quant_paths$quant_path, salmon_quant_paths$sampleID)
isoformExpression <- importIsoformExpression(
    sampleVector      = quant_vec,
    ignoreAfterPeriod = TRUE,   # match importRdata setting; handles versioned transcript IDs
    showProgress      = FALSE,
    quiet             = FALSE
)

# Step 2 — build SwitchAnalyzeRlist combining expression + annotation
# addAnnotatedORFs = FALSE: ORF extraction is handled explicitly in ISAR_EXTRACT_ORF
# to avoid running addORFfromGTF twice.
cat("Building SwitchAnalyzeRlist...\n")
switchAnalyzeRlist <- importRdata(
    isoformCountMatrix   = isoformExpression$counts,
    isoformRepExpression = isoformExpression$abundance,
    designMatrix         = design,
    isoformExonAnnoation = opt$gtf,
    isoformNtFasta       = opt$transcript_fasta,
    addAnnotatedORFs     = FALSE,
    ignoreAfterPeriod    = TRUE,   # Salmon IDs may carry transcript version suffixes
    showProgress         = FALSE
)

cat("\nImport summary:\n")
cat("  Genes:", length(unique(switchAnalyzeRlist$isoformFeatures$gene_id)), "\n")
cat("  Isoforms:", nrow(switchAnalyzeRlist$isoformFeatures), "\n")

# Save
saveRDS(switchAnalyzeRlist, file = opt$output)
cat("\nSaved:", opt$output, "\n")
cat("===========================================\n")
