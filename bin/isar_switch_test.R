#!/usr/bin/env Rscript

# Statistical testing for isoform switches using satuRn

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

# Apply defaults for optional numeric args
if (is.null(opt$alpha))      opt$alpha      <- "0.05"
if (is.null(opt$dif_cutoff)) opt$dif_cutoff <- "0.1"
opt$alpha      <- as.numeric(opt$alpha)
opt$dif_cutoff <- as.numeric(opt$dif_cutoff)

if (is.null(opt$input) || is.null(opt$output)) {
    stop("Required arguments: --input --output")
}

cat("===========================================\n")
cat("IsoformSwitchAnalyzeR - Switch Test\n")
cat("===========================================\n")

# Load data
switchAnalyzeRlist <- readRDS(opt$input)

cat("Loaded data:\n")
cat("  Genes:", length(unique(switchAnalyzeRlist$isoformFeatures$gene_id)), "\n")
cat("  Isoforms:", nrow(switchAnalyzeRlist$isoformFeatures), "\n")

# Pre-filter to remove low-expressed isoforms
cat("\nPre-filtering isoforms...\n")
switchAnalyzeRlist <- preFilter(
    switchAnalyzeRlist = switchAnalyzeRlist,
    geneExpressionCutoff = 1,       # min FPKM/TPM
    isoformExpressionCutoff = 0,    # no minimum for isoform
    IFcutoff = 0.01,                # min 1% isoform fraction
    removeSingleIsoformGenes = TRUE,
    reduceToSwitchingGenes = FALSE  # keep all, filter later
)

cat("After filtering:\n")
cat("  Genes:", length(unique(switchAnalyzeRlist$isoformFeatures$gene_id)), "\n")
cat("  Isoforms:", nrow(switchAnalyzeRlist$isoformFeatures), "\n")

# Test for isoform switches using satuRn
cat("\nRunning isoform switch test (satuRn)...\n")
cat("  Alpha:", opt$alpha, "\n")
cat("  dIF cutoff:", opt$dif_cutoff, "\n")

switchAnalyzeRlist <- tryCatch({
    isoformSwitchTestSatuRn(
        switchAnalyzeRlist = switchAnalyzeRlist,
        reduceToSwitchingGenes = TRUE,
        alpha = opt$alpha,
        dIFcutoff = opt$dif_cutoff
    )
}, error = function(e) {
    if (grepl("No genes were considered switching", conditionMessage(e))) {
        # No significant switches at these cutoffs — run without reduction
        # so test statistics are preserved in isoformFeatures
        cat("  NOTE: No significant switches at current thresholds, saving full results\n")
        isoformSwitchTestSatuRn(
            switchAnalyzeRlist = switchAnalyzeRlist,
            reduceToSwitchingGenes = FALSE,
            alpha = opt$alpha,
            dIFcutoff = opt$dif_cutoff
        )
    } else {
        stop(e)
    }
})

cat("\nSwitch test results:\n")
if (!is.null(switchAnalyzeRlist$isoformSwitchAnalysis) &&
    nrow(switchAnalyzeRlist$isoformSwitchAnalysis) > 0) {
    cat("  Switching isoforms:", nrow(switchAnalyzeRlist$isoformSwitchAnalysis), "\n")
    cat("  Switching genes:",
        length(unique(switchAnalyzeRlist$isoformSwitchAnalysis$gene_id)), "\n")
} else {
    cat("  No significant switches detected\n")
}

# Save
saveRDS(switchAnalyzeRlist, file = opt$output)
cat("\nSaved:", opt$output, "\n")
cat("===========================================\n")
