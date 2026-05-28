# LeafCutter container
# Python 3 + R with leafcutter bioconductor package
# Based on bioconductor docker for consistency with isar/report containers

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.description="LeafCutter for intron usage quantification and differential splicing"
LABEL org.opencontainers.image.source="https://github.com/BCI-KRP/alternative-splicing-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

# System tools: samtools + build deps for regtools
# regtools (Griffith Lab) is not in any standard apt repo — must be compiled from source
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    python3 \
    python3-pip \
    git \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && git clone --depth 1 https://github.com/griffithlab/regtools.git /opt/regtools \
    && cmake -S /opt/regtools -B /opt/regtools/build \
    && cmake --build /opt/regtools/build \
    && install -m755 /opt/regtools/build/regtools /usr/local/bin/ \
    && rm -rf /opt/regtools \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies for leafcutter scripts
RUN pip3 install --no-cache-dir --break-system-packages \
    numpy \
    scipy \
    pysam

# Pre-install ALL of leafcutter's DESCRIPTION dependencies before the GitHub build.
# remotes::install_github resolves CRAN deps but not Bioconductor ones (TailRank),
# so everything must be in place before the C++ compilation step.

# CRAN Imports + LinkingTo from leafcutter DESCRIPTION
RUN R -e " \
    install.packages(c( \
        'rstan', 'rstantools', 'StanHeaders', 'BH', 'Rcpp', 'RcppEigen', \
        'optparse', 'R.utils', 'ggplot2', 'reshape2', 'gridExtra', \
        'Hmisc', 'dplyr', 'foreach', 'magrittr', 'stringr', 'doMC', \
        'shiny', 'shinyjs', 'DT', 'intervals', 'shinycssloaders', 'gtable' \
    ), repos='https://cloud.r-project.org', Ncpus=4) \
"

# Bioconductor Imports: DirichletMultinomial + TailRank
# TailRank is Bioconductor-only; remotes cannot install it from CRAN repos
RUN R -e " \
    BiocManager::install(c('DirichletMultinomial', 'TailRank'), update=FALSE, ask=FALSE) \
"

# Install leafcutter from GitHub — all deps already satisfied, only C++ compilation needed
RUN R -e " \
    if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes') \
    remotes::install_github('davidaknowles/leafcutter/leafcutter', upgrade='never') \
"

# Verify leafcutter loads
RUN R -e "library(leafcutter); sessionInfo()"

WORKDIR /data
