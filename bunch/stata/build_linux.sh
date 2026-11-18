#!/bin/sh
set -eu
# Written but not executed here. Docker storage is caller-managed.
# Mount sources read-only; all package outputs, downloads, caches go under B.
repo=$(CDPATH= cd -- "$(dirname "$0")/../.." && pwd)
build="${SPECPKG_BUILD:-${TMPDIR:-/tmp}/spectral-build}/stata"
mkdir -p "$build/tmp" "$build/cache"
docker run --rm --platform linux/amd64 \
  -v "$repo:/src:ro" -v "$build:/build" \
  -e HOME=/build -e TMPDIR=/build/tmp -e XDG_CACHE_HOME=/build/cache \
  -e CC=gcc -e CXX=g++ \
  quay.io/pypa/manylinux_2_28_x86_64 \
  /opt/python/cp311-cp311/bin/python -B /src/bunch/stata/build.py \
  --platform linux64 --build-root /build
