# IsoformSwitchAnalyzeR Container
# Based on Bioconductor 3.20 (R 4.4.x)
# Includes all dependencies for isoform switch detection and reporting

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.description="IsoformSwitchAnalyzeR for isoform switch detection and functional annotation"
LABEL org.opencontainers.image.source="https://github.com/BCI-KRP/alternative-splicing-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

# Bioconductor packages: core ISAR stack
RUN R -e " \
    BiocManager::install(c( \
        'IsoformSwitchAnalyzeR', \
        'satuRn', \
        'DEXSeq', \
        'tximeta', \
        'tximport', \
        'BiocParallel' \
    ), update=FALSE, ask=FALSE) \
"

# CRAN packages used in R scripts and reports
RUN R -e " \
    install.packages(c( \
        'dplyr', 'tidyr', 'ggplot2', 'DT', \
        'rmarkdown', 'knitr', 'plotly', 'UpSetR' \
    ), repos='https://cloud.r-project.org', Ncpus=4) \
"

# pfamAnalyzeR — required for PFAM domain analysis (future)
RUN R -e "BiocManager::install('remotes', update=FALSE, ask=FALSE)" && \
    R -e "remotes::install_github('kvittingseerup/pfamAnalyzeR', upgrade='never')"

# Verify core packages load
RUN R -e "library(IsoformSwitchAnalyzeR); library(satuRn); library(tximeta); sessionInfo()"

WORKDIR /data
