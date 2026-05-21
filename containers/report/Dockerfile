# Report Container
# Based on rocker/verse which includes R, tidyverse, rmarkdown, and LaTeX
# Pre-installs all packages needed by render_report.Rmd — no runtime installs

FROM rocker/verse:4.4.2

LABEL org.opencontainers.image.description="R Markdown report renderer for alternative-splicing-nextflow"
LABEL org.opencontainers.image.source="https://github.com/BCI-KRP/alternative-splicing-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

# Interactive table and plot packages not bundled in rocker/verse
RUN R -e " \
    install.packages(c( \
        'DT', 'plotly', 'UpSetR' \
    ), repos='https://cloud.r-project.org', Ncpus=4) \
"

# Verify packages load
RUN R -e "library(DT); library(plotly); library(UpSetR); library(rmarkdown)"

WORKDIR /data
