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

# Clone leafcutter and regenerate Stan C++ model code for current rstan/StanHeaders.
# Leafcutter ships Stan 2.19.1 pre-compiled headers; rstan 2.21+ replaced the
# model base class (prob_grad -> model_base_crtp), causing compilation failure.
# rstan_config() regenerates src/stanExports_*.h, src/stanExports_*.cc, and
# R/stanmodels.R from the .stan source files using the installed rstan version.
RUN git clone --depth 1 https://github.com/davidaknowles/leafcutter.git /opt/leafcutter-src \
    && Rscript -e "setwd('/opt/leafcutter-src/leafcutter'); rstantools::rstan_config()"

RUN R CMD INSTALL /opt/leafcutter-src/leafcutter

# Verify the exact entrypoints used by the pipeline.
RUN python3 /opt/leafcutter-src/clustering/leafcutter_cluster_regtools.py -h >/dev/null \
    && ds_script="$(Rscript -e 'cat(system.file("scripts", "leafcutter_ds.R", package="leafcutter"))')" \
    && Rscript "${ds_script}" -h >/dev/null

WORKDIR /data
