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

# Add ORFs from GTF annotation; track success to gate analyzeNovelIsoformORF
cat("\nAdding ORFs from GTF...\n")
orf_loaded <- FALSE
switchAnalyzeRlist <- tryCatch({
    result <- addORFfromGTF(
        switchAnalyzeRlist = switchAnalyzeRlist,
        pathToGTF = opt$gtf,
        overwriteExistingORF = TRUE
    )
    orf_loaded <<- TRUE
    result
}, error = function(e) {
    if (grepl("No ORFs could be added", conditionMessage(e), fixed = TRUE)) {
        # Transcript IDs may carry version suffixes not present in the GTF — retry.
        cat("  NOTE: Retrying ORF mapping with ignoreAfterPeriod=TRUE\n")
        tryCatch({
            result <- addORFfromGTF(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToGTF = opt$gtf,
                overwriteExistingORF = TRUE,
                ignoreAfterBar = TRUE,
                ignoreAfterSpace = TRUE,
                ignoreAfterPeriod = TRUE
            )
            orf_loaded <<- TRUE
            result
        }, error = function(e2) {
            if (grepl("No ORFs could be added", conditionMessage(e2), fixed = TRUE)) {
                # GTF has no CDS entries (e.g. StringTie output, test data).
                # analyzeNovelIsoformORF v2.6 requires prior ORF annotation —
                # skip it and proceed directly to sequence export.
                cat("  NOTE: No CDS in GTF — skipping ORF analysis\n")
                switchAnalyzeRlist
            } else {
                stop(e2)
            }
        })
    } else {
        stop(e)
    }
})

# analyzeNovelIsoformORF requires base ORF annotation from addORFfromGTF (ISAR v2.6).
# Only run when ORFs were successfully loaded above.
if (orf_loaded) {
    cat("\nAnalyzing novel isoform ORFs...\n")
    switchAnalyzeRlist <- tryCatch({
        analyzeNovelIsoformORF(
            switchAnalyzeRlist = switchAnalyzeRlist,
            analysisAllIsoformsWithoutORF = TRUE
        )
    }, error = function(e) {
        cat("  NOTE: Skipping novel ORF inference:", conditionMessage(e), "\n")
        switchAnalyzeRlist
    })
} else {
    cat("  NOTE: Skipping analyzeNovelIsoformORF (no ORF annotation available)\n")
}

# Export sequences for external annotation tools
cat("\nExporting sequences...\n")
tryCatch({
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
}, error = function(e) {
    # Sequence export can fail when no ORFs/sequences are available (e.g. test data).
    # The RDS is saved regardless; FASTA outputs are declared optional in the module.
    cat("  NOTE: Sequence export failed (no FASTAs will be staged):", conditionMessage(e), "\n")
})

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
