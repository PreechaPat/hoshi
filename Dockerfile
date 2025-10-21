FROM ghcr.io/astral-sh/uv:python3.13-bookworm AS builder
WORKDIR /workspace

# Only keep fresh build artifacts inside the container.
RUN rm -rf dist

COPY pyproject.toml uv.lock README.md ./
COPY src src
COPY templates templates
COPY data data

RUN uv build --wheel

FROM python:3.13.7-slim-trixie

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /app

COPY --from=builder /workspace/dist/*.whl /tmp/dist/

RUN python -m pip install --no-cache-dir /tmp/dist/* \
    && rm -rf /tmp/dist

CMD ["hoshi"]
