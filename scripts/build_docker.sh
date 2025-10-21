#!/usr/bin/env bash
set -euo pipefail

# Resolve repository root relative to this script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

if ! command -v uv >/dev/null 2>&1; then
  echo "Error: uv is not installed or not on PATH." >&2
  echo "Install uv from https://github.com/astral-sh/uv and try again." >&2
  exit 1
fi

rm -rf dist/
uv build --wheel

if command -v python3 >/dev/null 2>&1; then
  PYTHON=python3
elif command -v python >/dev/null 2>&1; then
  PYTHON=python
else
  echo "Error: Neither python3 nor python is available on PATH." >&2
  exit 1
fi

# Read the package version from pyproject.toml.
PACKAGE_VERSION="$(
  "${PYTHON}" - <<'PY'
import tomllib
from pathlib import Path

data = tomllib.loads(Path("pyproject.toml").read_text())
print(data["project"]["version"])
PY
)"

if [[ -z "${PACKAGE_VERSION}" ]]; then
  echo "Error: Unable to determine package version from pyproject.toml." >&2
  exit 1
fi

VERSION_TAG="hoshi:${PACKAGE_VERSION}"
LATEST_TAG="hoshi:latest"
docker build -t "${VERSION_TAG}" -t "${LATEST_TAG}" .

echo "Built Docker image tags:"
echo "  ${VERSION_TAG}"
echo "  ${LATEST_TAG}"
