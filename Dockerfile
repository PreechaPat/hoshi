FROM ghcr.io/astral-sh/uv:python3.13-bookworm AS builder
WORKDIR /workspace

# Only keep fresh build artifacts inside the container.
RUN rm -rf dist

COPY pyproject.toml uv.lock README.md ./
COPY src src

RUN uv build --wheel

FROM python:3.13.7-slim-trixie

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Add ps, since nextflow need it.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps

WORKDIR /app

COPY --from=builder /workspace/dist/*.whl /tmp/dist/

RUN python -m pip install --no-cache-dir /tmp/dist/* \
    && rm -rf /tmp/dist

# Final cleanup
RUN apt-get clean \
 && rm -rf /var/lib/apt/lists/* /var/cache/apt/*

CMD ["hoshi"]
