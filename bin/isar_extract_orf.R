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
switchAnalyzeRlist <- tryCatch({
    addORFfromGTF(
        switchAnalyzeRlist = switchAnalyzeRlist,
        pathToGTF = opt$gtf,
        overwriteExistingORF = TRUE
    )
}, error = function(e) {
    if (grepl("No ORFs could be added", conditionMessage(e), fixed = TRUE)) {
        # Common mismatch: transcript IDs in quantification include version suffixes
        # while GTF matching requires IDs truncated before the period.
        cat("  NOTE: Retrying ORF mapping with ignoreAfterPeriod=TRUE\n")
        addORFfromGTF(
            switchAnalyzeRlist = switchAnalyzeRlist,
            pathToGTF = opt$gtf,
            overwriteExistingORF = TRUE,
            ignoreAfterBar = TRUE,
            ignoreAfterSpace = TRUE,
            ignoreAfterPeriod = TRUE
        )
    } else {
        stop(e)
    }
})

# Analyze novel isoform ORFs (those not in GTF)
cat("\nAnalyzing novel isoform ORFs...\n")
switchAnalyzeRlist <- tryCatch({
    analyzeNovelIsoformORF(
        switchAnalyzeRlist = switchAnalyzeRlist,
        analysisAllIsoformsWithoutORF = TRUE
    )
}, error = function(e) {
    if (grepl("genomeObject argument must be supplied", conditionMessage(e), fixed = TRUE)) {
        cat("  NOTE: Skipping novel ORF inference (missing genomeObject and/or transcript sequences in object)\n")
        switchAnalyzeRlist
    } else {
        stop(e)
    }
})

# Export sequences for external annotation tools
cat("\nExporting sequences...\n")
if (exists("exportSequences", mode = "function")) {
    exportSequences(
        switchAnalyzeRlist = switchAnalyzeRlist,
        pathToOutput = opt$output_dir,
        writeTranscriptSequences = TRUE,
        writeORFSequences = TRUE,
        writePeptideSequences = TRUE
    )
} else if (exists("extractSequence", mode = "function")) {
    switchAnalyzeRlist <- extractSequence(
        switchAnalyzeRlist = switchAnalyzeRlist,
        onlySwitchingGenes = FALSE,
        extractNTseq = TRUE,
        extractAAseq = TRUE,
        writeToFile = TRUE,
        pathToOutput = opt$output_dir,
        outputPrefix = "isoformSwitchAnalyzeR",
        quiet = FALSE
    )
} else {
    stop("No compatible sequence export function found (expected exportSequences or extractSequence)")
}

fasta_files <- list.files(opt$output_dir, pattern = "^isoformSwitchAnalyzeR_.*\\.fasta$", full.names = TRUE)
cat("Exported FASTA files:\n")
if (length(fasta_files) > 0) {
    for (f in fasta_files) {
        cat("  -", f, "\n")
    }
} else {
    cat("  - No FASTA files found matching isoformSwitchAnalyzeR_*.fasta\n")
}

# Save
saveRDS(switchAnalyzeRlist, file = opt$output)
cat("\nSaved:", opt$output, "\n")
cat("===========================================\n")
