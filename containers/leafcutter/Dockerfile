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

# R packages: leafcutter + dependencies
# leafcutter is installed from GitHub (davidaknowles/leafcutter)
RUN R -e " \
    install.packages(c( \
        'optparse', 'R.utils', 'rstan', 'ggplot2', 'reshape2', \
        'gridExtra', 'Hmisc', 'dplyr' \
    ), repos='https://cloud.r-project.org', Ncpus=4) \
"

RUN R -e " \
    BiocManager::install('DirichletMultinomial', update=FALSE, ask=FALSE) \
"

RUN R -e " \
    if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes') \
    remotes::install_github('davidaknowles/leafcutter/leafcutter', upgrade='never') \
"

# Verify leafcutter loads
RUN R -e "library(leafcutter); sessionInfo()"

WORKDIR /data
