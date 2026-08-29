#!/usr/bin/env Rscript
# Diagnose satuRn q-value pinning (uses the *_tested.rds outputs of ISAR_SWITCH_TEST).
#
# Usage:
#   Rscript diagnose_isar.R [path/to/*_tested.rds ...]
#   Rscript diagnose_isar.R <dir>            # globs **/*_tested.rds under dir
#   Rscript diagnose_isar.R                  # globs the default CRUK work dir
# Note: deliberately does NOT load the IsoformSwitchAnalyzeR package — the import
# only needs readRDS + base R, keeping the memory footprint small.

default_work <- "/gpfs/scratch/qp241615/altsplicing-cruk-batch_2026_04"
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  files <- list.files(default_work, pattern = "_tested\\.rds$",
                      recursive = TRUE, full.names = TRUE)
} else if (length(args) == 1 && dir.exists(args[1])) {
  files <- list.files(args[1], pattern = "_tested\\.rds$",
                      recursive = TRUE, full.names = TRUE)
} else {
  files <- args
}
files <- unique(files)

if (length(files) == 0) {
  stop("No *_tested.rds files found")
}

cat("Diagnosing", length(files), "ISAR_SWITCH_TEST RDS files\n\n")

for (f in files) {
  sar <- readRDS(f)
  # Access the switchAnalyzeRlist contents without loading IsoformSwitchAnalyzeR:
  # strip the S4 wrapper so the slots become a plain named list.
  sar <- unclass(sar)
  cat("=== ", basename(f), " ===\n", sep = "")
  cat("file:", f, "\n")

  feat <- sar$isoformFeatures
  if (!"isoform_switch_q_value" %in% colnames(feat)) {
    cat("  no isoform_switch_q_value column; skipping\n\n")
    next
  }
  q <- feat$isoform_switch_q_value[!is.na(feat$isoform_switch_q_value)]
  d <- feat$dIF[!is.na(feat$isoform_switch_q_value)]

  cat("n tested:", length(q), "\n")
  if (length(q) > 0) {
    cat("q range:", paste(range(q), collapse = " "),
        " IQR:", IQR(q), "\n")
  } else {
    cat("q range: NA IQR: NA\n")
  }
  cat("n |dIF| >= 0.1:", sum(abs(d) >= 0.1, na.rm = TRUE),
      " max |dIF|:", if (length(d) > 0) max(abs(d), na.rm = TRUE) else NA, "\n")
  cat("max |dIF| non-sig:",
      ifelse(any(q >= 0.05),
             max(abs(d[q >= 0.05]), na.rm = TRUE), NA), "\n")

  cat("designMatrix:\n")
  dm <- sar$designMatrix
  if (!is.null(dm)) {
    nr <- min(nrow(dm), 60)
    print(head(dm, nr))
    lv <- unique(dm$condition)
    cat("condition levels:", length(lv), "-", paste(lv, collapse = ", "), "\n")
    if (length(lv) == 2) {
      feat_lv <- unique(c(feat$condition_1, feat$condition_2))
      cat("condition values in isoformFeatures:",
          paste(feat_lv, collapse = ", "), "\n")
      cat("designMatrix levels == isoformFeatures levels:",
          identical(sort(as.character(lv)), sort(as.character(feat_lv))), "\n")
      # Sample order must match the expression matrix columns (minus the key
      # column), otherwise condition assignments are shifted.
      expr <- sar$isoformRepExpression
      if (!is.null(expr)) {
        expr_samples <- setdiff(colnames(expr), "isoform_id")
        cat("designMatrix sampleID order == expression columns order:",
            identical(as.character(dm$sampleID), as.character(expr_samples)), "\n")
      }
    }
  } else {
    cat("  (designMatrix is NULL)\n")
  }
  cat("\n")
}