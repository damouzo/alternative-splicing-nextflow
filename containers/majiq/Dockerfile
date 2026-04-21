# MAJIQ V3 Custom Container
# Academic license required - obtain from majiq.biociphers.org
# License file injected at runtime via MAJIQ_LICENSE_FILE env var

FROM python:3.12-slim

LABEL description="MAJIQ V3 for Local Splicing Variation analysis"
LABEL maintainer="BCI-KRP"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    libhts-dev \
    zlib1g-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install MAJIQ from Bitbucket
RUN pip install --no-cache-dir git+https://bitbucket.org/biociphers/majiq_academic.git

# Verify installation
RUN majiq --version || echo "MAJIQ installed (license required at runtime)"

WORKDIR /data

# Note: MAJIQ_LICENSE_FILE environment variable must be set at runtime
# pointing to a valid majiq_license*.key file
