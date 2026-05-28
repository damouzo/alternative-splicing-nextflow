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

# pfamAnalyzeR — required by analyzePFAM() for domain gain/loss detection
RUN R -e "BiocManager::install('remotes', update=FALSE, ask=FALSE)" && \
    R -e "remotes::install_github('kvittingseerup/pfamAnalyzeR', upgrade='never')"

# HMMER for PFAM domain analysis (analyzePFAM)
RUN apt-get update && apt-get install -y --no-install-recommends hmmer && \
    rm -rf /var/lib/apt/lists/*

# Python3 for scripting; iupred3 is not on PyPI (academic software, requires
# registration at iupred3.elte.hu) — analyzeIUPred2A() uses the web API by default
RUN apt-get update && apt-get install -y --no-install-recommends python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Verify core packages load
RUN R -e "library(IsoformSwitchAnalyzeR); library(satuRn); library(tximeta); sessionInfo()"

WORKDIR /data
