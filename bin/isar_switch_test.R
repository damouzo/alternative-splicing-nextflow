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
if (is.null(opt$alpha))            opt$alpha            <- "0.05"
if (is.null(opt$dif_cutoff))       opt$dif_cutoff       <- "0.1"
if (is.null(opt$gene_expr_cutoff)) opt$gene_expr_cutoff <- "1"
if (is.null(opt$iso_expr_cutoff))  opt$iso_expr_cutoff  <- "1"
if (is.null(opt$test_method))      opt$test_method      <- "auto"
opt$alpha            <- as.numeric(opt$alpha)
opt$dif_cutoff       <- as.numeric(opt$dif_cutoff)
opt$gene_expr_cutoff <- as.numeric(opt$gene_expr_cutoff)
opt$iso_expr_cutoff  <- as.numeric(opt$iso_expr_cutoff)
opt$test_method      <- tolower(opt$test_method)

if (!opt$test_method %in% c("auto", "dexseq", "satuRn")) {
    stop("Unknown --test_method '", opt$test_method,
         "' (expected 'auto', 'dexseq' or 'satuRn')")
}

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
    geneExpressionCutoff = opt$gene_expr_cutoff,       # min FPKM/TPM at gene level
    isoformExpressionCutoff = opt$iso_expr_cutoff,     # min TPM per isoform
    IFcutoff = 0.01,                                   # min 1% isoform fraction
    removeSingleIsoformGenes = TRUE,
    reduceToSwitchingGenes = FALSE  # keep all, filter later
)

cat("After filtering:\n")
cat("  Genes:", length(unique(switchAnalyzeRlist$isoformFeatures$gene_id)), "\n")
cat("  Isoforms:", nrow(switchAnalyzeRlist$isoformFeatures), "\n")

# Resolve the switch-test engine. 'auto' follows the adequacy rule used by
# ISAR's high-level workflow: the group with the fewest replicates limits the
# degrees of freedom, so DEXSeq (LRT, negative-binomial) is selected when the
# smallest condition has <= 5 replicates, and satuRn otherwise. Isolate the
# effective method here so it can be recorded verbatim in the outcome.
resolve_test_method <- function(sar, requested) {
    if (requested != "auto") return(requested)
    nr <- sar$conditions$nrReplicates
    nr <- nr[!is.na(nr)]
    if (length(nr) == 0) {
        warning("Could not resolve replicate counts for --test_method auto; defaulting to DEXSeq")
        return("dexseq")
    }
    if (min(nr) <= 5) "dexseq" else "satuRn"
}

test_method <- resolve_test_method(switchAnalyzeRlist, opt$test_method)
writeLines(test_method, con = "effective_test_method.txt")

# Test for isoform switches. DEXSeq (LRT, negative binomial) is preferred at low
# replication: SatuRn's empirical FDR (locfdr) step is known to collapse all
# q-values to one pinned value in high-isoform runs, so it is reserved (via
# --test_method satuRn or 'auto' with a well-replicated smallest group) for
# designs where it has enough power.
cat("\nRunning isoform switch test...\n")
cat("  Requested method:", opt$test_method, "\n")
cat("  Effective method:", test_method, "\n")
cat("  Alpha:", opt$alpha, "\n")
cat("  dIF cutoff:", opt$dif_cutoff, "\n")

run_switch_test <- function(sar, reduce_only) {
    if (test_method == "dexseq") {
        isoformSwitchTestDEXSeq(
            switchAnalyzeRlist      = sar,
            reduceToSwitchingGenes  = reduce_only,
            alpha                   = opt$alpha,
            dIFcutoff               = opt$dif_cutoff
        )
    } else {
        isoformSwitchTestSatuRn(
            switchAnalyzeRlist      = sar,
            reduceToSwitchingGenes  = reduce_only,
            alpha                   = opt$alpha,
            dIFcutoff               = opt$dif_cutoff
        )
    }
}

switchAnalyzeRlist <- tryCatch({
    run_switch_test(switchAnalyzeRlist, TRUE)
}, error = function(e) {
    if (grepl("No genes were considered switching", conditionMessage(e))) {
        # No significant switches at these cutoffs — run without reduction
        # so test statistics are preserved in isoformFeatures
        cat("  NOTE: No significant switches at current thresholds, saving full results\n")
        run_switch_test(switchAnalyzeRlist, FALSE)
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
