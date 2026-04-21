# IsoformSwitchAnalyzeR Container
# Based on Bioconductor 3.19 (R 4.3.x)

FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL description="IsoformSwitchAnalyzeR for isoform switch detection and functional annotation"
LABEL maintainer="BCI-KRP"

# Install IsoformSwitchAnalyzeR and dependencies
RUN R -e "BiocManager::install('IsoformSwitchAnalyzeR', update=FALSE, ask=FALSE)"

# Install pfamAnalyzeR (required for PFAM domain analysis)
RUN R -e "BiocManager::install('remotes', update=FALSE, ask=FALSE)" && \
    R -e "remotes::install_github('kvittingseerup/pfamAnalyzeR', upgrade='never')"

# Install additional dependencies for functional annotation
RUN R -e "BiocManager::install(c('satuRn', 'DEXSeq', 'tximeta'), update=FALSE, ask=FALSE)"

# Verify installation
RUN R -e "library(IsoformSwitchAnalyzeR); sessionInfo()"

WORKDIR /data
