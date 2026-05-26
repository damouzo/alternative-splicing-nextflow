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

consequences_to_analyze <- c(
    'intron_retention',
    'coding_potential',
    'ORF_seq_similarity',
    'NMD_status'
)

# Ensure intron retention has been precomputed when requested by ISAR.
switchAnalyzeRlist <- tryCatch({
    analyzeIntronRetention(
        switchAnalyzeRlist = switchAnalyzeRlist,
        onlySwitchingGenes = FALSE,
        dIFcutoff = opt$dif_cutoff,
        showProgress = FALSE,
        quiet = TRUE
    )
}, error = function(e) {
    cat("  NOTE: Skipping intron_retention analysis:", conditionMessage(e), "\n")
    consequences_to_analyze <<- consequences_to_analyze[consequences_to_analyze != 'intron_retention']
    switchAnalyzeRlist
})

# Analyze consequences with graceful degradation when optional annotations are missing
no_switching_genes <- FALSE
while (length(consequences_to_analyze) > 0) {
    analysis_attempt <- tryCatch({
        list(
            ok = TRUE,
            result = analyzeSwitchConsequences(
                switchAnalyzeRlist = switchAnalyzeRlist,
                consequencesToAnalyze = consequences_to_analyze,
                dIFcutoff = opt$dif_cutoff,
                onlySigIsoforms = FALSE,
                showProgress = FALSE
            )
        )
    }, error = function(e) {
        list(ok = FALSE, message = conditionMessage(e))
    })

    if (analysis_attempt$ok) {
        switchAnalyzeRlist <- analysis_attempt$result
        break
    }

    err_msg <- analysis_attempt$message
    if (grepl("coding_potential", err_msg, fixed = TRUE) &&
        'coding_potential' %in% consequences_to_analyze) {
        cat("  NOTE: Skipping coding_potential (missing CPAT/CPC2 results)\n")
        consequences_to_analyze <- consequences_to_analyze[consequences_to_analyze != 'coding_potential']
    } else if (grepl("No genes were considered switching", err_msg, fixed = TRUE)) {
        cat("  NOTE: No genes considered switching with current cutoffs; skipping consequence analysis\n")
        no_switching_genes <- TRUE
        break
    } else if (grepl("intron retention", tolower(err_msg), fixed = TRUE) &&
               'intron_retention' %in% consequences_to_analyze) {
        cat("  NOTE: Skipping intron_retention (classification not available)\n")
        consequences_to_analyze <- consequences_to_analyze[consequences_to_analyze != 'intron_retention']
    } else {
        stop(err_msg)
    }
}

if (length(consequences_to_analyze) == 0) {
    cat("  NOTE: No consequence categories available with current annotations\n")
}

# Extract top switches
cat("\nExtracting top switches...\n")
if (!is.null(switchAnalyzeRlist$isoformSwitchAnalysis) &&
    nrow(switchAnalyzeRlist$isoformSwitchAnalysis) > 0) {
    topSwitches <- tryCatch({
        extractTopSwitches(
            switchAnalyzeRlist = switchAnalyzeRlist,
            filterForConsequences = FALSE,  # Include all switches
            n = Inf,
            extractGenes = FALSE,
            sortByQvals = TRUE
        )
    }, error = function(e) {
        cat("  NOTE: No significant switching isoforms available for top switch extraction\n")
        data.frame()
    })
    
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
consequenceSummary <- tryCatch({
    if (no_switching_genes) {
        stop("No switching genes available for consequence summary")
    }
    extractConsequenceSummary(
        switchAnalyzeRlist = switchAnalyzeRlist,
        includeCombined = TRUE,
        consequencesToPlot = 'all',
        asFractionTotal = FALSE
    )
}, error = function(e) {
    cat("  NOTE: Could not build consequence summary:", conditionMessage(e), "\n")
    data.frame()
})

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
