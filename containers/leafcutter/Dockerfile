# LeafCutter container
# Python 3 + R with leafcutter bioconductor package
# Based on bioconductor docker for consistency with isar/report containers

FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.description="LeafCutter for intron usage quantification and differential splicing"
LABEL org.opencontainers.image.source="https://github.com/damouzo/alternative-splicing-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

# System tools: samtools + build deps for regtools
# regtools (Griffith Lab) is not in any standard apt repo — must be compiled from source
RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    python3 \
    python3-pip \
    git \
    cmake \
    libtbb-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    && if [ ! -f /usr/include/tbb/tbb_stddef.h ] && [ -f /usr/include/oneapi/tbb/version.h ]; then \
        mkdir -p /usr/include/tbb; \
        printf '#pragma once\n#include <oneapi/tbb/version.h>\n' > /usr/include/tbb/tbb_stddef.h; \
    fi \
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

# Pin the Stan toolchain. rstan delegates Stan parsing to the stanc.js bundled
# in StanHeaders, so a newer StanHeaders (>= 2.33, current CRAN) rejects
# leafcutter's pre-2.33 array syntax (`int ys[N]`) and rstan_config() aborts.
# Versions fixed to the image that produced the published results.
RUN R -e " \
    if (!requireNamespace('remotes', quietly=TRUE)) install.packages('remotes', repos='https://cloud.r-project.org', Ncpus=4); \
    remotes::install_version('StanHeaders', version='2.32.10', repos='https://cloud.r-project.org', upgrade='never', Ncpus=4); \
    remotes::install_version('rstan',       version='2.32.7',  repos='https://cloud.r-project.org', upgrade='never', Ncpus=4); \
    remotes::install_version('rstantools',  version='2.6.0',   repos='https://cloud.r-project.org', upgrade='never', Ncpus=4) \
"

# Remaining CRAN Imports + LinkingTo from leafcutter DESCRIPTION.
# rstan/StanHeaders/rstantools are pinned above and must not appear here.
RUN R -e " \
    install.packages(c( \
        'BH', 'Rcpp', 'RcppEigen', \
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

# Clone leafcutter at the exact commit used to build the validated image and
# regenerate Stan C++ model code with the pinned toolchain. Pinning the commit
# keeps .stan sources and model code identical to the published results.
RUN git clone https://github.com/davidaknowles/leafcutter.git /opt/leafcutter-src \
    && git -C /opt/leafcutter-src checkout 2c9907ef66adf0bfb3092f0ceb6886ee5c046fbb \
    && Rscript -e "setwd('/opt/leafcutter-src/leafcutter'); rstantools::rstan_config()"

RUN R CMD INSTALL /opt/leafcutter-src/leafcutter

# Verify the exact entrypoints used by the pipeline.
# leafcutter_ds.R lives in the repo's scripts/ dir, not in the R package inst/scripts/,
# so system.file() returns "" — use the source path directly.
RUN python3 /opt/leafcutter-src/clustering/leafcutter_cluster_regtools.py -h >/dev/null \
    && Rscript /opt/leafcutter-src/scripts/leafcutter_ds.R --help >/dev/null

WORKDIR /data
