FROM ghcr.io/astral-sh/uv:0.11.8-python3.13-trixie AS builder

WORKDIR /workspace

# Only keep fresh build artifacts inside the container.
RUN rm -rf dist

COPY pyproject.toml uv.lock README.md ./
COPY src src

RUN uv build --wheel

# production
FROM python:3.13.7-slim-trixie

ARG VERSION
LABEL org.opencontainers.image.version=${VERSION}
LABEL maintainer="Preecha Patumcharoenpol"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install runtime dependencies:
# - procps: provides `ps`, required by nextflow
# - libpango/libgdk-pixbuf/libffi/libharfbuzz: native libraries required by
#   WeasyPrint for PDF generation (libgobject/libglib come in transitively)
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
    procps \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libgdk-pixbuf-2.0-0 \
    libffi8 \
    libharfbuzz0b \
    fonts-dejavu-core

WORKDIR /app

COPY --from=builder /workspace/dist/*.whl /tmp/dist/

RUN python -m pip install --no-cache-dir /tmp/dist/* \
    && rm -rf /tmp/dist

# Final cleanup
RUN apt-get clean \
 && rm -rf /var/lib/apt/lists/* /var/cache/apt/*

CMD ["hoshi"]
