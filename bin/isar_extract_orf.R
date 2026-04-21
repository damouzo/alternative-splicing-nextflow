#!/usr/bin/env Rscript

# Extract ORFs from transcript sequences

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

if (is.null(opt$output_dir)) opt$output_dir <- "."

if (is.null(opt$input) || is.null(opt$gtf) || is.null(opt$output)) {
    stop("Required arguments: --input --gtf --output")
}

cat("===========================================\n")
cat("IsoformSwitchAnalyzeR - Extract ORF\n")
cat("===========================================\n")

# Load data
switchAnalyzeRlist <- readRDS(opt$input)

cat("Loaded data with", nrow(switchAnalyzeRlist$isoformFeatures), "isoforms\n")

# Add ORFs from GTF annotation
cat("\nAdding ORFs from GTF...\n")
switchAnalyzeRlist <- addORFfromGTF(
    switchAnalyzeRlist = switchAnalyzeRlist,
    pathToGTF = opt$gtf,
    overwriteExistingORF = TRUE
)

# Analyze novel isoform ORFs (those not in GTF)
cat("\nAnalyzing novel isoform ORFs...\n")
switchAnalyzeRlist <- analyzeNovelIsoformORF(
    switchAnalyzeRlist = switchAnalyzeRlist,
    analysisAllIsoformsWithoutORF = TRUE
)

# Export sequences for external annotation tools
cat("\nExporting sequences...\n")
exportSequences(
    switchAnalyzeRlist = switchAnalyzeRlist,
    pathToOutput = opt$output_dir,
    writeTranscriptSequences = TRUE,
    writeORFSequences = TRUE,
    writePeptideSequences = TRUE
)

cat("Exported FASTA files:\n")
cat("  -", file.path(opt$output_dir, "isoformSwitchAnalyzeR_isoform.fasta"), "\n")
cat("  -", file.path(opt$output_dir, "isoformSwitchAnalyzeR_ORF.fasta"), "\n")
cat("  -", file.path(opt$output_dir, "isoformSwitchAnalyzeR_AA.fasta"), "\n")

# Save
saveRDS(switchAnalyzeRlist, file = opt$output)
cat("\nSaved:", opt$output, "\n")
cat("===========================================\n")
