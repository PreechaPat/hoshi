#!/usr/bin/env bash
set -euo pipefail

# Get the short 10-character git commit hash
SHA_SHORT=$(git rev-parse --short=10 HEAD)
SHA_FULL=$(git rev-parse HEAD)

# Uses the layer cache. Pass --no-cache (or --pull) to force fresh apt
# packages / base image, e.g.: ./build-rev.sh --no-cache
echo "Building hoshi image with revision ${SHA_SHORT}"
docker build . -t "hoshi:sha-${SHA_SHORT}" -t "hoshi:latest" --build-arg "VERSION=${SHA_SHORT}" "$@"
