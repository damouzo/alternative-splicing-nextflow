FROM python:3.11-slim

LABEL maintainer="BCI-KRP pipeline" \
      description="PEGASAS (Python 3 port) + R dependencies"

# System libraries required by R and Python scientific stack
RUN apt-get update && apt-get install -y --no-install-recommends \
        r-base \
        r-base-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        build-essential \
        pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies
RUN pip install --no-cache-dir numpy scipy matplotlib

# Install PEGASAS Python 3 package
WORKDIR /opt/pegasas
COPY . .
RUN pip install --no-cache-dir .

# R packages: LSD and data.table (required by cor_matrix_direct_perm.R)
RUN Rscript -e "install.packages(c('LSD', 'data.table', 'ggplot2'), repos='https://cloud.r-project.org', Ncpus=4)"

# Verify installation
RUN python3 -c "import PEGASAS; print('PEGASAS import OK')" && \
    PEGASAS --version

WORKDIR /data
