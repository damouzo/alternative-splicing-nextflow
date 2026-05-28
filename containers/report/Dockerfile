# Report Container
# Based on rocker/verse which includes R, tidyverse, rmarkdown, and LaTeX
# Pre-installs all packages needed by render_report.Rmd — no runtime installs

FROM rocker/verse:4.4.2

LABEL org.opencontainers.image.description="R Markdown report renderer for alternative-splicing-nextflow"
LABEL org.opencontainers.image.source="https://github.com/BCI-KRP/alternative-splicing-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

# Interactive table and visualization packages not bundled in rocker/verse
RUN R -e " \
    install.packages(c( \
        'DT', 'plotly', 'UpSetR', 'ggrepel' \
    ), repos='https://cloud.r-project.org', Ncpus=4) \
"

# BiocManager is required to install Bioconductor packages
RUN R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Bioconductor packages for GO/KEGG enrichment and pathway analysis
# org.Hs.eg.db / org.Mm.eg.db: annotation databases for human and mouse
# clusterProfiler: GO and KEGG over-representation and gene-set enrichment
# fgsea: fast GSEA (used for PEGASAS-equivalent per-sample pathway scoring)
# update=FALSE removed: rocker/verse ships pre-installed CRAN packages that may
# be outdated relative to what clusterProfiler requires; allowing updates ensures
# the full dependency graph resolves cleanly
RUN R -e " \
    BiocManager::install( \
        c('clusterProfiler', 'org.Hs.eg.db', 'org.Mm.eg.db', 'fgsea'), \
        ask=FALSE \
    ) \
"

# Verify all packages load cleanly
RUN R -e " \
    library(DT); library(plotly); library(UpSetR); library(ggrepel); \
    library(rmarkdown); library(clusterProfiler); library(fgsea) \
"

WORKDIR /data
