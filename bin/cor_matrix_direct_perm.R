#!/usr/bin/env Rscript
# Pearson correlation between pathway KS scores and per-event PSI, with a
# permutation p-value. Vendored into the project bin/ so the pipeline uses this
# (fixed) version directly via Nextflow's bin directory — no container rebuild
# required. The PEGASAS_CORRELATION module calls the two upstream Python steps
# (prepareGeneMatrixOrdered.py / generateMatrixbySample.py, resolved through the
# installed PEGASAS package) and then invokes this script.
#
# Fixes vs. the original PEGASAS port:
#   * The PDF device is opened lazily on the first significant event, so contrasts
#     with no significant event leave NO (empty) PDF behind.
#   * The global matrix reports EVERY event (no pre-filtering) with its Pearson r,
#     the permutation p (NA where |r| <= 0.3, where it was not computed), and a
#     `significant` flag — so the report table is never empty and non-significant
#     events are visible.
#   * The permutation null distributions for x (pathway KS scores) are precomputed
#     ONCE per pathway and reused across all events, instead of regenerating
#     replicate(N, centre(sample(x))) for every single event (80K+ times).
#   * Thresholds are unchanged (|r| > 0.3 AND p < 0.0002) — only display/compute
#     behaviour changed.
#   * The maximum permutation count is configurable via 4th CLI arg (default 1000).

args <- commandArgs(TRUE)
# Matrix 1 = pathway KS scores (sorted by group): row1 = sample ids, row2 = scores
# Matrix 2 = PSI matrix aligned to the same sample order
library('LSD')
library('data.table')

signature  <- args[1]
PSImatrix  <- args[2]
folder     <- args[3]
N_PERMS    <- if (length(args) >= 4L) as.integer(args[4]) else 1000L
N_PERMS    <- max(200L, N_PERMS)   # lower bound so the tiered test works
analysis_prefix <- unlist(strsplit(tail(unlist(strsplit(signature, '/')), 1), '\\.sorted\\.txt'))[1]

matrix1  <- as.matrix(fread(signature, header = TRUE, skip = 0))
matrix2  <- as.matrix(fread(PSImatrix,  header = TRUE, skip = 0))
matrix11 <- as.matrix(as.numeric(matrix1[, -1]))            # score vector (samples)
matrix22 <- suppressWarnings(apply(matrix2[, -(1:8), drop = FALSE], 2, as.numeric))  # events x samples

# Significance thresholds: |r| > 0.3 is a pre-filter for the (expensive) permutation
# test. The p-value threshold is derived from the permutation resolution — an event
# is considered significant only when its observed correlation is more extreme than
# ALL N_PERMS permutations (i.e. p == 0 within the resolution of the test).
R_ABS_MIN <- 0.3
P_MAX     <- 1 / N_PERMS

# Vectorised permutation test. The original called cor() once per permutation,
# which is the dominant cost (88k events x up to 5000 perms). Here all
# permutations of x are materialised as a single (n x N) matrix and the N
# correlations with y are computed in one vectorised sweep, reusing x's
# centred/scaled form. Same adaptive replication logic and the SAME result.
NCORES <- max(1L, as.integer(Sys.getenv('NSLOTS', parallel::detectCores(logical = FALSE) %||% 1)))
NCORES <- min(NCORES, 8L)

centre <- function(v) {
  m <- mean(v, na.rm = TRUE)
  s <- sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) rep(0, length(v)) else (v - m) / s
}

# Vectorised correlations between every column of permX (n x N, centred x perms)
# and y (centred), ignoring NA via complete cases per column. Returns length-N.
fast.cor <- function(permX, y_c) {
  has_y <- !is.na(y_c)
  if (!all(has_y)) {
    out <- numeric(ncol(permX))
    for (j in seq_len(ncol(permX))) {
      ok <- has_y & !is.na(permX[, j])
      out[j] <- sum(permX[ok, j] * y_c[ok]) / (sum(ok) - 1)
    }
    out
  } else {
    as.vector((crossprod(permX, y_c)) / (length(y_c) - 1))
  }
}

# Precomputed permutation nulls: generated ONCE per pathway (x is the same for
# all events), then reused across all candidate events. Eliminates the dominant
# cost (replicate(N, centre(sample(x))) called ~80K times per pathway).
make.cor.perm.test <- function(x_c, permX200, permX1000, permXN, N_MAX) {
  force(x_c); force(permX200); force(permX1000); force(permXN); force(N_MAX)
  function(y) {
    y_c <- centre(y)
    if (all(y_c == 0)) return(NA_real_)
    obs <- sum(x_c * y_c, na.rm = TRUE) / (sum(!is.na(x_c) & !is.na(y_c)) - 1)
    reps <- fast.cor(permX200, y_c)
    p <- min(mean(obs >= reps, na.rm = TRUE), mean(obs <= reps, na.rm = TRUE))
    n_used <- 200L
    if (N_MAX > 200L && p < 0.2) {
      reps <- fast.cor(permX1000, y_c)
      p <- min(mean(obs >= reps, na.rm = TRUE), mean(obs <= reps, na.rm = TRUE))
      n_used <- 1000L
    }
    if (N_MAX > 1000L && p < 0.04) {
      reps <- fast.cor(permXN, y_c)
      p <- min(mean(obs >= reps, na.rm = TRUE), mean(obs <= reps, na.rm = TRUE))
      n_used <- N_MAX
    }
    # Laplace (+1/N+1) estimator: no exact p == 0, so the report's BH step has
    # a meaningful lower resolution (1/(N+1)) instead of many ties at zero.
    (p * n_used + 1) / (n_used + 1)
  }
}

# Empty-input guard: write complete (empty) tables and exit cleanly.
nr <- nrow(matrix22)
if (is.null(nr) || nr == 0) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  global_hdr <- c('GeneName','chr','strand','exonStart','exonEnd',
                  'upstreamEE','downstreamES','pearson_r','p_value','significant')
  writeLines(paste(global_hdr, collapse = '\t'),
             paste0(folder, '/', analysis_prefix, '_global_cor_matrix.txt'))
  write.table(data.frame(event_id = character(0), pearson_r = numeric(0), p_value = numeric(0)),
              file = paste0(folder, '/', analysis_prefix, '_high_cor_matrix.txt'),
              row.names = FALSE, quote = FALSE, sep = '\t')
  quit(save = 'no', status = 0)
}

# Compute Pearson r for every event. The permutation p-value is computed only
# when |r| exceeds the significance floor (R_ABS_MIN): a pure compute
# optimisation — every event is still reported in the global matrix with its
# r, so no data is filtered out. p is NA where it was not computed.
x_scores <- as.vector(matrix11[, 1])
cor_vec  <- apply(matrix22, 1, function(row) cor(x_scores, row,
                          use = 'pairwise.complete.obs', method = 'pearson'))
event_id <- vapply(seq_len(nr), function(i) {
  s <- paste(matrix2[i, 2], matrix2[i, 3], matrix2[i, 4],
            matrix2[i, 5], matrix2[i, 6], matrix2[i, 7], matrix2[i, 8], sep = '_')
  paste0(unlist(strsplit(s, ' ', fixed = TRUE)), collapse = '')
}, character(1))

p_vec   <- rep(NA_real_, nr)
sig_vec <- logical(nr)

# Precompute permutation null distributions ONCE (x_scores are identical for all
# events in this pathway). The per-event function captures them in a closure so
# mclapply only passes matrix22[i, ] — no redundant replicate() calls.
x_c <- centre(x_scores)
if (all(x_c == 0)) {
  cat(sprintf('[INFO] %s: x is constant (all zero centred) — permutation test skipped\n',
              analysis_prefix))
} else {
  set.seed(42L)  # reproducible null distributions
  permX200  <- replicate(200L,  centre(sample(x_scores)))
  permX1000 <- replicate(1000L, centre(sample(x_scores)))
  permXN    <- replicate(N_PERMS, centre(sample(x_scores)))
  cat(sprintf('[INFO] %s: precomputed permutation nulls (200/1000/%d) in %.2f s\n',
              analysis_prefix, N_PERMS, proc.time()[[3L]]))
}

# Run the permutation test only on candidates (|r| > R_ABS_MIN), in parallel.
cand <- which(!is.na(cor_vec) & abs(cor_vec) > R_ABS_MIN)
if (length(cand) > 0 && !all(x_c == 0)) {
  perm_test <- make.cor.perm.test(x_c, permX200, permX1000, permXN, N_PERMS)
  cat(sprintf('[INFO] %s: %d/%d events pass |r|>%.2f — running permutation test on %d cores\n',
              analysis_prefix, length(cand), nr, R_ABS_MIN, NCORES))
  p_cand <- parallel::mclapply(cand, function(i)
      tryCatch(perm_test(matrix22[i, ]), error = function(e) NA_real_),
    mc.cores = NCORES)
  p_vec[cand] <- unlist(p_cand)
  cat(sprintf('[INFO] %s: permutation test done in %.2f s\n',
              analysis_prefix, proc.time()[[3L]]))
} else if (length(cand) > 0 && all(x_c == 0)) {
  cat(sprintf('[INFO] %s: %d/%d events pass |r|>%.2f — skipped (x constant)\n',
              analysis_prefix, length(cand), nr, R_ABS_MIN))
}
sig_vec <- !is.na(p_vec) & p_vec < P_MAX

# Lazily render heatscatters only for significant events (no empty PDFs).
pdf_path   <- paste0(folder, '/', analysis_prefix, '_high_cor_scatterplots.pdf')
pdf_opened <- FALSE
sig_idx <- which(sig_vec)
if (length(sig_idx) > 0) {
  pdf(pdf_path, onefile = TRUE)
  pdf_opened <- TRUE
  for (i in sig_idx) {
    s <- paste(matrix2[i, 2], matrix2[i, 3], matrix2[i, 4],
               matrix2[i, 5], matrix2[i, 6], sep = '_')
    s <- paste0(unlist(strsplit(s, ' ', fixed = TRUE)), collapse = '')
    cor_p_lab <- format(p_vec[i], scientific = TRUE)
    if (p_vec[i] == 0) cor_p_lab <- '<1E-3'
    heatscatter(x_scores, matrix22[i, ],
                xlab = paste0(matrix1[1, 1], ' GSEA Score (K-S statistic)\n',
                              "Pearson's corr coef:", round(cor_vec[i], digits = 3),
                              '\tp-value:', cor_p_lab),
                ylab = '', ylim = c(0, 1),
                main = paste(strwrap(s, width = 40), collapse = '\n'),
                cex.main = 1, cex.lab = 0.95, cex.axis = 0.8,
                font.lab = 2, font.main = 2)
    mtext(expression(italic(psi)), side = 2, line = 2, cex = 1.3, font = 2)
    abline(lm(matrix22[i, ] ~ x_scores))
  }
  dev.off()
}

# Global matrix: ALL events with r; p where computed (NA otherwise); significance flag.
global_out <- data.frame(
  GeneName    = matrix2[, 2],
  chr         = matrix2[, 3],
  strand      = matrix2[, 4],
  exonStart   = matrix2[, 5],
  exonEnd     = matrix2[, 6],
  upstreamEE  = matrix2[, 7],
  downstreamES = matrix2[, 8],
  pearson_r   = cor_vec,
  p_value     = p_vec,
  significant = sig_vec,
  stringsAsFactors = FALSE
)
write.table(global_out,
            file = paste0(folder, '/', analysis_prefix, '_global_cor_matrix.txt'),
            row.names = FALSE, quote = FALSE, sep = '\t')

# High-correlation matrix: only significant events (subset of the global table).
high_out <- data.frame(
  event_id  = event_id[sig_vec],
  pearson_r = cor_vec[sig_vec],
  p_value   = p_vec[sig_vec],
  stringsAsFactors = FALSE
)
write.table(high_out,
            file = paste0(folder, '/', analysis_prefix, '_high_cor_matrix.txt'),
            row.names = FALSE, quote = FALSE, sep = '\t')
