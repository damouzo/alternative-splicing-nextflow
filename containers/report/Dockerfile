# Report Container
# Based on bioconductor_docker:RELEASE_3_20 instead of rocker/verse to avoid
# dependency conflicts between rocker/verse's CRAN snapshot and clusterProfiler
# from Bioc 3.20. bioconductor_docker has BiocManager 3.20 pre-configured and no
# conflicting pre-installed CRAN packages.

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.description="R Markdown report renderer for alternative-splicing-nextflow"
LABEL org.opencontainers.image.source="https://github.com/BCI-KRP/alternative-splicing-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

# pandoc is required by rmarkdown for HTML rendering
RUN apt-get update && apt-get install -y --no-install-recommends pandoc \
    && rm -rf /var/lib/apt/lists/*

# R Markdown stack + visualization packages
RUN R -e " \
    install.packages(c( \
        'rmarkdown', 'knitr', 'ggplot2', 'dplyr', \
        'DT', 'plotly', 'UpSetR', 'ggrepel' \
    ), repos='https://cloud.r-project.org', Ncpus=4) \
"

# Bioconductor packages for GO/KEGG enrichment and pathway analysis
# org.Hs.eg.db / org.Mm.eg.db: annotation databases for human and mouse
# clusterProfiler: GO and KEGG over-representation and gene-set enrichment
# fgsea: fast GSEA (used for PEGASAS-equivalent per-sample pathway scoring)
RUN R -e " \
    BiocManager::install( \
        c('clusterProfiler', 'org.Hs.eg.db', 'org.Mm.eg.db', 'fgsea'), \
        ask=FALSE, update=FALSE \
    ) \
"

# Verify all packages load cleanly
RUN R -e " \
    library(DT); library(plotly); library(UpSetR); library(ggrepel); \
    library(rmarkdown); library(clusterProfiler); library(fgsea) \
"

WORKDIR /data
