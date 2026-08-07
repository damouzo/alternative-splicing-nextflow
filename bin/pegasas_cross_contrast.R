#!/usr/bin/env Rscript
# Cross-contrast PEGASAS comparison: which pathways harbour significantly
# correlated splicing events, and how those pathways overlap across contrasts.
#
# Consumes one <comp_id>_sig_pathways.tsv per contrast (pathway, n_sig_events),
# produced by the PEGASAS_COLLATE step, and writes:
#   upset_significant_pathways.png  — UpSet plot of significant-pathway sets
#   heatmap_pathway_contrast.png    — pathway x contrast count heatmap
#   cross_contrast_summary.tsv       — pathway x contrast (n_sig_events) long table
#
# Contrasts with zero significant events still appear (columns of zeros), so the
# comparison is never silently empty.

suppressPackageStartupMessages({
  library(ggplot2)
  library(UpSetR)
})

args <- commandArgs(TRUE)
out_dir <- if (length(args) >= 1L && nzchar(args[1])) args[1] else '.'
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

draw_text_png <- function(path, lines) {
  png(path, width = 800, height = 400)
  plot.new()
  text(0.5, 0.5, paste(lines, collapse = '\n'), cex = 1.1)
  dev.off()
}

sig_files <- list.files('.', pattern = '_sig_pathways\\.tsv$', full.names = TRUE)
contrasts <- sub('_sig_pathways\\.tsv$', '', basename(sig_files))

if (length(sig_files) == 0L) {
  draw_text_png(file.path(out_dir, 'upset_significant_pathways.png'),
                'No contrasts found for PEGASAS cross-contrast comparison.')
  draw_text_png(file.path(out_dir, 'heatmap_pathway_contrast.png'),
                'No contrasts found for PEGASAS cross-contrast comparison.')
  writeLines('pathway\tcontrast\tn_sig_events',
             file.path(out_dir, 'cross_contrast_summary.tsv'))
  quit(save = 'no', status = 0)
}

# Parse every per-contrast file into a long table (pathway, contrast, n_sig_events).
long_df <- do.call(rbind, lapply(sig_files, function(f) {
  comp_id <- sub('_sig_pathways\\.tsv$', '', basename(f))
  x <- tryCatch(read.table(f, header = TRUE, sep = '\t', stringsAsFactors = FALSE),
                error = function(e) NULL)
  if (is.null(x) || nrow(x) == 0L ||
      !all(c('pathway', 'n_sig_events') %in% names(x))) {
    # Keep a 0-row frame so do.call(rbind, ...) stays well-formed: every
    # column must be length 0 (a length-1 scalar would error in data.frame).
    return(data.frame(pathway = character(0), contrast = character(0),
                      n_sig_events = integer(0), stringsAsFactors = FALSE))
  }
  data.frame(pathway = x$pathway, contrast = comp_id,
             n_sig_events = as.integer(x$n_sig_events), stringsAsFactors = FALSE)
}))

all_pathways <- sort(unique(long_df$pathway))

# Wide matrix pathway x contrast (0 where a pathway is absent for a contrast).
build_wide <- function() {
  if (length(all_pathways) == 0L) return(NULL)
  m <- matrix(0L, nrow = length(all_pathways), ncol = length(contrasts),
             dimnames = list(all_pathways, contrasts))
  if (nrow(long_df) > 0L) {
    for (k in seq_len(nrow(long_df))) {
      m[long_df$pathway[k], long_df$contrast[k]] <- long_df$n_sig_events[k]
    }
  }
  m
}
wide <- build_wide()

# Long summary (always written, even when empty).
if (!is.null(wide)) {
  summary_long <- do.call(rbind, lapply(contrasts, function(cc) {
    data.frame(pathway = rownames(wide), contrast = cc,
               n_sig_events = wide[, cc], stringsAsFactors = FALSE)
  }))
} else {
  summary_long <- data.frame(pathway = character(0), contrast = character(0),
                             n_sig_events = integer(0), stringsAsFactors = FALSE)
}
write.table(summary_long, file.path(out_dir, 'cross_contrast_summary.tsv'),
            row.names = FALSE, quote = FALSE, sep = '\t')

# --- UpSet plot of significant-pathway sets across contrasts ---
sets <- lapply(contrasts, function(cc) {
  if (is.null(wide)) character(0) else rownames(wide)[wide[, cc] > 0L]
})
names(sets) <- contrasts

png(file.path(out_dir, 'upset_significant_pathways.png'),
    width = max(800, 60 * length(contrasts) + 400), height = 600)
if (any(lengths(sets) > 0L)) {
  upset(fromList(sets), order.by = 'freq',
        mainbar.y.label = 'Pathway intersections',
        sets.x.label = 'Significant pathways per contrast',
        text.scale = 1.1)
} else {
  plot.new()
  text(0.5, 0.5,
       'No pathway reached the significant-correlation threshold in any contrast.')
}
dev.off()

# --- Heatmap pathway x contrast (n_sig_events) ---
if (is.null(wide) || nrow(wide) == 0L) {
  draw_text_png(file.path(out_dir, 'heatmap_pathway_contrast.png'),
                'No pathways to display (no significant events).')
} else {
  totals <- rowSums(wide)
  pathway_order <- names(sort(totals, decreasing = TRUE))
  heat_df <- summary_long
  heat_df$pathway  <- factor(heat_df$pathway, levels = pathway_order)
  heat_df$contrast <- factor(heat_df$contrast, levels = contrasts)
  p <- ggplot(heat_df, aes(x = contrast, y = pathway, fill = n_sig_events)) +
    geom_tile(colour = 'white') +
    scale_fill_gradient(low = '#f7fbff', high = '#08306b', name = 'Sig. events') +
    labs(title = 'Pathway-correlated splicing events across contrasts',
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          panel.grid = element_blank())
  n_path <- length(pathway_order)
  ggsave(file.path(out_dir, 'heatmap_pathway_contrast.png'), p,
         width = max(6, 0.25 * length(contrasts) + 4),
         height = max(4, 0.22 * n_path + 2), dpi = 130, limitsize = FALSE)
}
