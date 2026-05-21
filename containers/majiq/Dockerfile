# MAJIQ V3 Container
# MAJIQ is distributed under an academic licence (https://majiq.biociphers.org).
# The licence file is NOT baked into this image — it must be provided at runtime
# via the MAJIQ_LICENSE_FILE environment variable.
#
# Build instructions:
#   1. Obtain your licence key from https://majiq.biociphers.org
#   2. docker build -t your-registry/majiq:3.0 containers/majiq/
#   3. docker push your-registry/majiq:3.0
#   4. Pass --majiq_sif your-registry/majiq:3.0 to the pipeline (Docker)
#      OR build the Apptainer .sif and pass --majiq_sif /path/to/majiq.sif

FROM python:3.12-slim

LABEL org.opencontainers.image.description="MAJIQ V3 for Local Splicing Variation analysis"
LABEL org.opencontainers.image.source="https://github.com/BCI-KRP/alternative-splicing-nextflow"

# Build dependencies for MAJIQ's C++/pybind11 extensions
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    libhts-dev \
    zlib1g-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Limit parallel compilation to avoid OOM during pybind11/CMake build
ENV CMAKE_BUILD_PARALLEL_LEVEL=2
ENV MAX_JOBS=2

# Install MAJIQ from the academic Bitbucket repository
RUN pip install --no-cache-dir git+https://bitbucket.org/biociphers/majiq_academic.git

# Verify the binary is available (will print version; licence checked at runtime)
RUN majiq-v3 --version || echo "MAJIQ installed — licence key required at runtime"

WORKDIR /data

# Runtime: set MAJIQ_LICENSE_FILE to the path of your majiq_license*.key file
# Example: docker run -e MAJIQ_LICENSE_FILE=/keys/majiq_license.key ...
