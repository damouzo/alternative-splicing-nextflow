# Sashimi plotting container
# Mirrors upstream rmats2sashimiplot environment and is published via GHCR.

FROM docker.io/xinglab/rmats2sashimiplot:latest

LABEL org.opencontainers.image.description="rmats2sashimiplot runtime for alternative-splicing-nextflow"
LABEL org.opencontainers.image.source="https://github.com/BCI-KRP/alternative-splicing-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

# rmats2sashimiplot calls samtools internally via subprocess.
# Install it explicitly and keep image slim.
RUN set -eux; \
    apt-get update; \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends samtools ca-certificates; \
    rm -rf /var/lib/apt/lists/*

# Fail fast at build time if required tools are unavailable.
RUN set -eux; \
    TOOL_BIN="$(command -v rmats2sashimiplot || true)"; \
    if [ -z "${TOOL_BIN}" ] && [ -x /usr/local/bin/rmats2sashimiplot ]; then \
        TOOL_BIN=/usr/local/bin/rmats2sashimiplot; \
    fi; \
    [ -n "${TOOL_BIN}" ]; \
    "${TOOL_BIN}" -h >/dev/null 2>&1 || "${TOOL_BIN}" --help >/dev/null 2>&1; \
    command -v samtools >/dev/null 2>&1; \
    samtools --version | head -1

WORKDIR /data
