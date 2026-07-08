# Sashimi plotting container
# Mirrors upstream rmats2sashimiplot environment and is published via GHCR.

FROM docker.io/xinglab/rmats2sashimiplot:latest

LABEL org.opencontainers.image.description="rmats2sashimiplot runtime for alternative-splicing-nextflow"
LABEL org.opencontainers.image.source="https://github.com/BCI-KRP/alternative-splicing-nextflow"
LABEL org.opencontainers.image.licenses="MIT"

# Fail fast at build time if rmats2sashimiplot is unavailable.
RUN set -eux; \
    TOOL_BIN="$(command -v rmats2sashimiplot || true)"; \
    if [ -z "${TOOL_BIN}" ] && [ -x /usr/local/bin/rmats2sashimiplot ]; then \
        TOOL_BIN=/usr/local/bin/rmats2sashimiplot; \
    fi; \
    [ -n "${TOOL_BIN}" ]; \
    "${TOOL_BIN}" -h >/dev/null 2>&1 || "${TOOL_BIN}" --help >/dev/null 2>&1

WORKDIR /data
