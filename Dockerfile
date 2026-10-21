FROM ghcr.io/astral-sh/uv:0.13.4-trixie-slim AS builder

WORKDIR /workspace

# Pinned runtime deps from the lockfile only, so this layer (and the pip
# install that uses it) is rebuilt only when uv.lock changes, not on every
# source edit.
COPY pyproject.toml uv.lock ./
RUN uv export --frozen --no-dev --no-emit-project --no-hashes -o requirements.txt

COPY README.md ./
COPY src src
RUN uv build --wheel

# production
FROM python:3.13.7-slim-trixie

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install runtime dependencies:
# - procps: provides `ps`, required by nextflow
# - libpango/libgdk-pixbuf/libffi/libharfbuzz: native libraries required by
#   WeasyPrint for PDF generation (libgobject/libglib come in transitively)
# Cleanup lives in the same RUN so the apt lists never land in a layer.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
    procps \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    libgdk-pixbuf-2.0-0 \
    libffi8 \
    libharfbuzz0b \
    fonts-dejavu-core \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/* /var/cache/apt/*

WORKDIR /app

COPY --from=builder /workspace/requirements.txt /tmp/requirements.txt
RUN python -m pip install --no-cache-dir -r /tmp/requirements.txt \
 && rm /tmp/requirements.txt

COPY --from=builder /workspace/dist/*.whl /tmp/dist/
RUN python -m pip install --no-cache-dir --no-deps /tmp/dist/*.whl \
 && rm -rf /tmp/dist

# Last, so a new VERSION (e.g. a git SHA per build) does not invalidate the
# cached layers above.
ARG VERSION
LABEL org.opencontainers.image.version=${VERSION}
LABEL maintainer="Preecha Patumcharoenpol"

CMD ["hoshi"]
