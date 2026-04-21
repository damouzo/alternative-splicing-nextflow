#!/usr/bin/env Rscript

# Analyze functional consequences of isoform switches

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

# Apply defaults for optional args
if (is.null(opt$output_dir))  opt$output_dir  <- "."
if (is.null(opt$dif_cutoff))  opt$dif_cutoff  <- "0.1"
if (is.null(opt$top_n_plots)) opt$top_n_plots <- "25"
opt$dif_cutoff  <- as.numeric(opt$dif_cutoff)
opt$top_n_plots <- as.integer(opt$top_n_plots)

if (is.null(opt$input) || is.null(opt$output)) {
    stop("Required arguments: --input --output")
}

cat("===========================================\n")
cat("IsoformSwitchAnalyzeR - Switch Consequences\n")
cat("===========================================\n")

# Load data
switchAnalyzeRlist <- readRDS(opt$input)

cat("Loaded data with", nrow(switchAnalyzeRlist$isoformFeatures), "isoforms\n")

# Note: External annotation tools (CPC2, PFAM, SignalP, etc.) would be run here
# For now, we analyze consequences based on available annotations (ORF-based)

cat("\nAnalyzing switch consequences...\n")
cat("  dIF cutoff:", opt$dif_cutoff, "\n")

# Analyze consequences
switchAnalyzeRlist <- analyzeSwitchConsequences(
    switchAnalyzeRlist = switchAnalyzeRlist,
    consequencesToAnalyze = c(
        'intron_retention',
        'coding_potential',
        'ORF_seq_similarity',
        'NMD_status'
        # Domain and signal peptide analysis require external tool results:
        # 'domains_identified', 'domain_isotype',
        # 'IDR_identified', 'IDR_type', 'IDR_seq_similarity',
        # 'signal_peptide_identified',
        # 'topology_identified'
    ),
    dIFcutoff = opt$dif_cutoff,
    onlySigIsoforms = FALSE,
    showProgress = FALSE
)

# Extract top switches
cat("\nExtracting top switches...\n")
if (!is.null(switchAnalyzeRlist$isoformSwitchAnalysis)) {
    topSwitches <- extractTopSwitches(
        switchAnalyzeRlist = switchAnalyzeRlist,
        filterForConsequences = FALSE,  # Include all switches
        n = Inf,
        extractGenes = FALSE,
        sortByQvals = TRUE
    )
    
    output_csv <- file.path(opt$output_dir, "top_isoform_switches.csv")
    write.csv(topSwitches, output_csv, row.names = FALSE)
    cat("  Saved:", output_csv, "\n")
    cat("  Total switches:", nrow(topSwitches), "\n")
} else {
    cat("  No switches detected\n")
    # Create empty file
    write.csv(data.frame(), file.path(opt$output_dir, "top_isoform_switches.csv"))
}

# Extract consequence summary
cat("\nExtracting consequence summary...\n")
consequenceSummary <- extractConsequenceSummary(
    switchAnalyzeRlist = switchAnalyzeRlist,
    includeCombined = TRUE,
    consequencesToPlot = 'all',
    asFractionTotal = FALSE
)

output_summary <- file.path(opt$output_dir, "consequence_summary.csv")
write.csv(consequenceSummary, output_summary, row.names = FALSE)
cat("  Saved:", output_summary, "\n")

# Generate switch plots for top N genes
if (!is.null(switchAnalyzeRlist$isoformSwitchAnalysis) && 
    nrow(switchAnalyzeRlist$isoformSwitchAnalysis) > 0) {
    
    cat("\nGenerating switch plots for top", opt$top_n_plots, "genes...\n")
    
    switchplot_dir <- file.path(opt$output_dir, "switchplots")
    dir.create(switchplot_dir, showWarnings = FALSE, recursive = TRUE)
    
    tryCatch({
        switchPlotTopSwitches(
            switchAnalyzeRlist = switchAnalyzeRlist,
            n = opt$top_n_plots,
            filterForConsequences = FALSE,
            fileType = "pdf",
            pathToOutput = switchplot_dir
        )
        cat("  Saved plots in:", switchplot_dir, "\n")
    }, error = function(e) {
        cat("  Warning: Could not generate switch plots:", e$message, "\n")
    })
}

# Save final object
saveRDS(switchAnalyzeRlist, file = opt$output)
cat("\nSaved final object:", opt$output, "\n")
cat("===========================================\n")
