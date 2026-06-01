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

# Install leafcutter from source and patch stanmodels.R for modern rstan API.
# Upstream still calls rstan::stanc_builder(), removed in recent rstan releases.
# Keep a fallback to rstan::stanc() so namespace loading works on current R images.
RUN git clone --depth 1 https://github.com/davidaknowles/leafcutter.git /opt/leafcutter-src \
    && Rscript - <<'RS'
f <- '/opt/leafcutter-src/leafcutter/R/stanmodels.R'
x <- readLines(f, warn = FALSE)
pat <- grep('stanfit <- rstan::stanc_builder\\(stan_file,', x)
if (length(pat) != 1) {
    stop('Expected stanc_builder call not found in stanmodels.R')
}
i <- pat[[1]]
x[i:(i+2)] <- c(
  '  stanc_fun <- if (exists("stanc_builder", where = asNamespace("rstan"), mode = "function")) {',
  '    rstan::stanc_builder',
  '  } else {',
  '    function(file, allow_undefined = TRUE, obfuscate_model_name = FALSE) {',
  '      rstan::stanc(file = file, allow_undefined = allow_undefined)',
  '    }',
  '  }',
  '  stanfit <- stanc_fun(stan_file,',
  '                       allow_undefined = TRUE,',
  '                       obfuscate_model_name = FALSE)'
)
writeLines(x, f)
RS
RUN R CMD INSTALL /opt/leafcutter-src/leafcutter

# Verify the exact entrypoints used by the pipeline.
RUN python3 /opt/leafcutter-src/clustering/leafcutter_cluster_regtools.py -h >/dev/null \
    && ds_script="$(Rscript -e 'cat(system.file("scripts", "leafcutter_ds.R", package="leafcutter"))')" \
    && Rscript "${ds_script}" -h >/dev/null

WORKDIR /data
