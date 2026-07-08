# Sashimi plotting container
# Mirrors upstream rmats2sashimiplot environment and is published via GHCR.

FROM docker.io/xinglab/rmats2sashimiplot:latest

LABEL org.opencontainers.image.description="rmats2sashimiplot runtime for alternative-splicing-nextflow"
LABEL org.opencontainers.image.source="https://github.com/BCI-KRP/alternative-splicing-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

# Fail fast at build time if entrypoint binary is missing.
RUN python3 /usr/local/bin/rmats2sashimiplot --version

WORKDIR /data
